# Fitted spatial covariance comparison with conditional Gaussian calibration.
# All geometry is expression-independent and reusable across pathways. The
# calibration preserves the complete score Gram, including FA correlations.

.passage_gp_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > .Machine$integer.max || x != floor(x))
    stop(name, " must be an integer >= ", lower)
  as.integer(x)
}

.passage_gp_matern <- function(a, b, range) {
  d2 <- outer(rowSums(a * a), rowSums(b * b), "+") - 2 * tcrossprod(a, b)
  z <- sqrt(3 * pmax(d2, 0)) / range
  (1 + z) * exp(-z)
}

# Deterministic farthest-point landmarks use coordinates alone.
.passage_gp_landmarks <- function(coords, m) {
  centered <- sweep(coords, 2L, colMeans(coords), "-")
  pick <- which.min(rowSums(centered^2))
  nearest <- rep(Inf, nrow(coords))
  selected <- integer()
  for (j in seq_len(min(m, nrow(coords)))) {
    selected <- c(selected, pick)
    delta <- sweep(coords, 2L, coords[pick, ], "-")
    nearest <- pmin(nearest, rowSums(delta^2))
    if (max(nearest) <= .Machine$double.eps * max(1, max(rowSums(centered^2)))) break
    pick <- which.max(nearest)
  }
  selected
}

passage_prepare_gp_kernel <- function(coords, X = NULL,
                                      ranges = c(0.05, 0.1, 0.2, 0.4),
                                      kernel_method = c("lowrank", "dense"),
                                      n_features = 64L, max_dense_n = 1000L) {
  started <- proc.time()[[3L]]
  kernel_method <- match.arg(kernel_method)
  coords <- passage_check_coords(coords)
  n <- nrow(coords)
  if (n < 5L) stop("at least five coordinates are required")
  ids <- rownames(coords)
  if (is.null(ids)) ids <- paste0("spot_", seq_len(n))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("coordinate IDs must be unique and nonempty")
  rownames(coords) <- ids
  if (!is.null(X) && !is.null(rownames(X)) && !identical(rownames(X), ids))
    stop("X spot IDs must exactly match coordinates in order")
  X <- passage_prepare_design(X, n)
  rownames(X) <- ids
  qx <- qr(X)
  df <- n - qx$rank
  if (qx$rank != ncol(X) || df < 3L) stop("X must have full rank and at least three residual degrees of freedom")
  if (!is.numeric(ranges) || !length(ranges) || any(!is.finite(ranges)) || any(ranges <= 0) || anyDuplicated(ranges))
    stop("ranges must be distinct positive finite numbers")
  ranges <- sort(as.numeric(ranges))
  n_features <- .passage_gp_integer(n_features, "n_features", 2L)
  max_dense_n <- .passage_gp_integer(max_dense_n, "max_dense_n", 5L)
  if (kernel_method == "dense" && n > max_dense_n)
    stop("dense GP kernel exceeds max_dense_n; use kernel_method='lowrank'")
  if (sum(apply(coords, 2L, stats::var)) <= 0) stop("coordinates have no spatial variation")
  landmarks <- if (kernel_method == "lowrank") .passage_gp_landmarks(coords, n_features) else integer()
  spectra <- vectors <- vector("list", length(ranges))
  for (j in seq_along(ranges)) {
    if (kernel_method == "dense") {
      R <- .passage_gp_matern(coords, coords, ranges[j])
      R <- qr.resid(qx, t(qr.resid(qx, R)))
      eig <- eigen((R + t(R)) / 2, symmetric = TRUE)
      positive <- eig$values > max(eig$values) * 1e-10
      lambda <- eig$values[positive]
      U <- eig$vectors[, positive, drop = FALSE]
    } else {
      W <- .passage_gp_matern(coords[landmarks, , drop = FALSE], coords[landmarks, , drop = FALSE], ranges[j])
      ew <- eigen((W + t(W)) / 2, symmetric = TRUE)
      positive <- ew$values > max(ew$values) * 1e-10
      C <- .passage_gp_matern(coords, coords[landmarks, , drop = FALSE], ranges[j])
      features <- C %*% sweep(ew$vectors[, positive, drop = FALSE], 2L, sqrt(ew$values[positive]), "/")
      features <- qr.resid(qx, features)
      sf <- svd(features, nu = min(dim(features)), nv = 0L)
      positive <- sf$d > max(sf$d) * 1e-8
      lambda <- sf$d[positive]^2
      U <- sf$u[, positive, drop = FALSE]
    }
    if (!length(lambda) || !all(is.finite(lambda)) || sum(lambda) <= .Machine$double.eps)
      stop("kernel has no residual spatial variation at range ", ranges[j])
    if (length(lambda) > df) stop("numerical kernel rank exceeds residual degrees of freedom")
    # Mean residual-space spatial variance is one: h is an average residual
    # variance fraction, including for the explicitly approximate GP.
    spectra[[j]] <- lambda * (df / sum(lambda))
    vectors[[j]] <- U
  }
  all_vectors <- do.call(cbind, vectors)
  # A single union basis lets null calibration operate in at most L*m
  # dimensions, irrespective of the number of spots. LAPACK QR is rank
  # revealing; eigen/SVD here also handles overlapping kernel eigenspaces.
  su <- svd(all_vectors, nu = min(ncol(all_vectors), df), nv = 0L)
  rank_union <- min(df, sum(su$d > max(su$d) * 1e-9))
  basis <- su$u[, seq_len(rank_union), drop = FALSE]
  projection_transform <- crossprod(basis, all_vectors)
  cursor <- 0L
  kernels <- lapply(seq_along(ranges), function(j) {
    ix <- cursor + seq_along(spectra[[j]])
    cursor <<- max(ix)
    list(range = ranges[j], lambda = spectra[[j]], indices = ix,
      transform = projection_transform[, ix, drop = FALSE])
  })
  ans <- list(coords = coords, X = X, df = df, basis = basis,
    projection_transform = projection_transform, kernels = kernels,
    ranges = ranges, kernel_method = kernel_method, n_features = n_features,
    landmarks = landmarks, normalization = "residual kernel trace equals residual df",
    preparation_seconds = proc.time()[[3L]] - started, cache_version = 1L)
  class(ans) <- c("passage_gp_kernel", "list")
  ans
}

# Profile over total scale analytically and over h/range on finite grids.
# Columns may be components from many null draws, evaluated in one BLAS batch.
.passage_gp_profile <- function(projection_sq, total_ss, cache, fractions) {
  projection_sq <- as.matrix(projection_sq)
  nc <- ncol(projection_sq)
  df <- cache$df
  best <- rep(0, nc)
  best_range <- rep(NA_real_, nc)
  best_h <- rep(0, nc)
  best_scale <- total_ss / df
  for (j in seq_along(cache$kernels)) {
    spec <- cache$kernels[[j]]
    energy <- projection_sq[spec$indices, , drop = FALSE]
    rest <- pmax(0, total_ss - colSums(energy))
    hs <- fractions[fractions > 0]
    if (!length(hs)) next
    v <- outer(spec$lambda - 1, hs, "*") + 1
    inverse_energy <- crossprod(1 / v, energy) + outer(1 / (1 - hs), rest)
    logdet <- colSums(log(v)) + (df - length(spec$lambda)) * log1p(-hs)
    gain <- df * (matrix(log(total_ss), nrow = length(hs), ncol = nc, byrow = TRUE) - log(inverse_energy)) - logdet
    choice <- max.col(t(gain), ties.method = "first")
    loc <- cbind(choice, seq_len(nc))
    value <- gain[loc]
    replace <- value > best
    best[replace] <- value[replace]
    best_range[replace] <- spec$range
    best_h[replace] <- hs[choice[replace]]
    best_scale[replace] <- inverse_energy[loc][replace] / df
  }
  fits <- data.frame(statistic = pmax(0, best), range = best_range,
    fraction = best_h, total_variance = best_scale,
    spatial_variance = best_h * best_scale,
    noise_variance = (1 - best_h) * best_scale,
    at_fraction_upper = best_h > 0 & best_h == max(fractions),
    at_range_boundary = !is.na(best_range) & best_range %in% range(cache$ranges))
  list(statistic = fits$statistic, fits = fits)
}

# First q rows of a Haar Stiefel frame, multiplied by chol(F'F).
# Bartlett's Wishart draw integrates out the d-q unused directions. The
# Cholesky has positive diagonal, avoiding the sign bias of qr.Q alone.
.passage_gp_rotation <- function(q, df, gram_chol) {
  K <- ncol(gram_chol)
  A <- matrix(stats::rnorm(q * K), q, K)
  rest_df <- df - q
  if (rest_df >= K) {
    extra <- stats::rWishart(1L, rest_df, diag(K))[, , 1L]
  } else if (rest_df > 0L) {
    z <- matrix(stats::rnorm(rest_df * K), rest_df, K)
    extra <- crossprod(z)
  } else extra <- matrix(0, K, K)
  upper <- chol(crossprod(A) + extra)
  # A %*% solve(upper) %*% gram_chol, without a general matrix inverse.
  A %*% backsolve(upper, gram_chol)
}

passage_test_pathway_gp <- function(fit, coords = NULL, kernel_cache = NULL,
                                    ranges = c(0.05, 0.1, 0.2, 0.4),
                                    fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
                                    B = 999L, seed = 1L,
                                    kernel_method = c("lowrank", "dense"),
                                    n_features = 64L, max_dense_n = 1000L,
                                    batch_size = 64L, keep_null = TRUE) {
  started <- proc.time()[[3L]]
  if (!inherits(fit, "passage_pathway_fit") ||
      !identical(fit$estimation_scope, "pathway_only") ||
      !identical(fit$scores_source, "unsmoothed_expression") ||
      !fit$method %in% c("pca", "factor_analysis"))
    stop("fit must be an unsmoothed passage_fit_pathway() result")
  scores <- as.matrix(fit$scores)
  if (!is.numeric(scores) || any(!is.finite(scores)) || ncol(scores) != fit$K)
    stop("fit has invalid scores")
  B <- .passage_gp_integer(B, "B")
  seed <- .passage_gp_integer(seed, "seed", 0L)
  batch_size <- .passage_gp_integer(batch_size, "batch_size")
  if (!is.logical(keep_null) || length(keep_null) != 1L || is.na(keep_null)) stop("keep_null must be TRUE or FALSE")
  if (!is.numeric(fractions) || !length(fractions) || any(!is.finite(fractions)) ||
      any(fractions < 0 | fractions >= 1) || !any(fractions == 0) || anyDuplicated(fractions))
    stop("fractions must be distinct numbers in [0,1), including 0")
  fractions <- sort(as.numeric(fractions))
  if (is.null(kernel_cache)) {
    if (is.null(coords)) stop("supply coords or a prepared kernel_cache")
    if (is.null(rownames(coords))) rownames(coords) <- rownames(scores)
    kernel_cache <- passage_prepare_gp_kernel(coords, X = fit$X, ranges = ranges,
      kernel_method = match.arg(kernel_method), n_features = n_features, max_dense_n = max_dense_n)
  } else {
    if (!inherits(kernel_cache, "passage_gp_kernel") || !identical(kernel_cache$cache_version, 1L))
      stop("invalid kernel_cache")
    if (!is.null(coords) && !isTRUE(all.equal(as.matrix(coords), kernel_cache$coords)))
      stop("coords conflict with kernel_cache")
    if (!missing(ranges) && !identical(sort(as.numeric(ranges)), kernel_cache$ranges))
      stop("ranges conflict with kernel_cache")
    if (!missing(kernel_method) && !identical(match.arg(kernel_method), kernel_cache$kernel_method))
      stop("kernel_method conflicts with kernel_cache")
    if (!missing(n_features) && n_features != kernel_cache$n_features)
      stop("n_features conflicts with kernel_cache")
  }
  cache <- kernel_cache
  if (!identical(rownames(scores), rownames(cache$coords)) ||
      !isTRUE(all.equal(fit$X, cache$X, check.attributes = TRUE)))
    stop("kernel_cache coordinates/design must exactly match the fitted pathway")
  K <- ncol(scores)
  df <- cache$df
  if (K < 1L || K >= df || qr(scores)$rank != K) stop("scores need full column rank below residual df")
  scores <- qr.resid(qr(cache$X), scores)
  total <- colSums(scores^2)
  if (any(total <= 0)) stop("zero residual component variance")
  # Unit-norm scores improve conditioning without changing any component LR.
  standardized <- sweep(scores, 2L, sqrt(total), "/")
  gram_chol <- tryCatch(chol(crossprod(standardized)), error = function(e)
    stop("score Gram is numerically singular"))
  projection <- crossprod(cache$projection_transform, crossprod(cache$basis, standardized))
  observed <- .passage_gp_profile(projection^2, rep(1, K), cache, fractions)
  statistic <- sum(observed$statistic)
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  null <- numeric(B)
  calibration_started <- proc.time()[[3L]]
  for (start in seq.int(1L, B, by = batch_size)) {
    nb <- min(batch_size, B - start + 1L)
    orientations <- matrix(0, ncol(cache$basis), K * nb)
    for (b in seq_len(nb)) orientations[, (b - 1L) * K + seq_len(K)] <-
      .passage_gp_rotation(ncol(cache$basis), df, gram_chol)
    pp <- crossprod(cache$projection_transform, orientations)
    values <- .passage_gp_profile(pp^2, rep(1, K * nb), cache, fractions)$statistic
    null[start + seq_len(nb) - 1L] <- colSums(matrix(values, K, nb))
  }
  # A numerical tolerance counts near-ties conservatively.
  exceedances <- sum(null >= statistic - 1e-10 * max(1, abs(statistic)))
  component_fits <- observed$fits
  for (nm in c("total_variance", "spatial_variance", "noise_variance"))
    component_fits[[nm]] <- component_fits[[nm]] * total
  component_fits <- cbind(component = colnames(scores), component_fits, stringsAsFactors = FALSE)
  list(statistic = statistic, p = (1 + exceedances) / (B + 1),
    component_fits = component_fits, null = if (keep_null) null else NULL,
    B = B, seed = seed, K = K, residual_df = df, exceedances = exceedances,
    mc_floor = 1 / (B + 1), method = "gp_profile_rotation",
    calibration = "conditional Gaussian rotation; complete score Gram preserved; range/fraction reselected in every null draw",
    null_assumption = "iid Gaussian residual rows conditional on fixed design, with arbitrary gene covariance; covariance-only pathway PCA/FA selection",
    role = "self-contained test of spatial covariance in the retained pathway components; not competitive",
    estimation = list(method = fit$method, scope = fit$estimation_scope, K = K,
      rank_rule = fit$rank_rule, gene_names = fit$gene_names),
    configuration = list(ranges = cache$ranges, fractions = fractions,
      kernel_method = cache$kernel_method, n_features_requested = cache$n_features,
      landmarks_used = length(cache$landmarks), kernel_ranks = lengths(lapply(cache$kernels, `[[`, "lambda")),
      union_rank = ncol(cache$basis), normalization = cache$normalization,
      batch_size = batch_size, coordinates = "supplied units; no automatic rescaling",
      statistic = "sum of per-component twice-profile-restricted-loglikelihood improvements; not a full joint FA likelihood"),
    timing = c(total_seconds = proc.time()[[3L]] - started,
      calibration_seconds = proc.time()[[3L]] - calibration_started,
      kernel_preparation_seconds = cache$preparation_seconds),
    interpretation = "finite-grid plug-in spatial fractions/ranges; uncertainty and power need study-specific assessment; lowrank uses an explicitly approximated covariance")
}
