# SpaPath: minimal subset-adaptive spatial pathway testing prototype.
#
# The first-paper target is deliberately narrow:
#   H0_any(G): no member of pathway G has residual spatial variability.
# It computes feature-level spatial variance-component score statistics,
# combines them with ASSET-like softmax/top-k pathway statistics, reports
# driver features, and estimates a simple pathway spatial variance fraction.

spapath_test <- function(Y,
                         coords,
                         pathways,
                         X = NULL,
                         ranges = NULL,
                         m = 20,
                         kernel = c("matern32", "gaussian", "exponential"),
                         alpha_grid = c(0.5, 1, 2, 4, 8),
                         topk_grid = NULL,
                         n_sim = 2000,
                         min_pathway_size = 2,
                         max_pathway_size = 500,
                         fdr_method = "BH",
                         seed = NULL,
                         verbose = TRUE) {
  kernel <- match.arg(kernel)
  Y <- .spapath_check_y(Y)
  coords <- .spapath_check_coords(coords, nrow(Y))
  X <- .spapath_prepare_design(X, nrow(Y))
  pathways <- .spapath_check_pathways(pathways, colnames(Y))

  if (is.null(ranges)) {
    ranges <- spapath_default_ranges(coords)
  }
  ranges <- as.numeric(ranges)
  if (any(!is.finite(ranges)) || any(ranges <= 0)) {
    stop("ranges must be positive finite numbers")
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (verbose) {
    message("SpaPath: computing feature-level spatial scores")
  }

  scores <- spapath_feature_scores(
    Y = Y,
    coords = coords,
    X = X,
    ranges = ranges,
    m = m,
    kernel = kernel
  )

  if (verbose) {
    message("SpaPath: testing ", length(pathways), " pathways")
  }

  pathway_results <- spapath_pathway_tests(
    scores = scores,
    pathways = pathways,
    alpha_grid = alpha_grid,
    topk_grid = topk_grid,
    n_sim = n_sim,
    min_pathway_size = min_pathway_size,
    max_pathway_size = max_pathway_size
  )

  pathway_results$fdr <- stats::p.adjust(pathway_results$p_value, method = fdr_method)
  pathway_results <- pathway_results[order(pathway_results$p_value), , drop = FALSE]
  rownames(pathway_results) <- NULL

  out <- list(
    results = pathway_results,
    feature_scores = scores$feature_scores,
    ranges = ranges,
    m = m,
    kernel = kernel,
    alpha_grid = alpha_grid,
    call = match.call()
  )
  class(out) <- "spapath_result"
  out
}

spapath_feature_scores <- function(Y,
                                   coords,
                                   X = NULL,
                                   ranges = NULL,
                                   m = 20,
                                   kernel = c("matern32", "gaussian", "exponential")) {
  kernel <- match.arg(kernel)
  Y <- .spapath_check_y(Y)
  coords <- .spapath_check_coords(coords, nrow(Y))
  X <- .spapath_prepare_design(X, nrow(Y))
  if (is.null(ranges)) {
    ranges <- spapath_default_ranges(coords)
  }

  fit0 <- .spapath_residualize(Y, X)
  resid <- fit0$resid
  sigma2 <- fit0$sigma2
  E <- sweep(resid, 2, sqrt(sigma2), "/")

  nn <- spapath_knn(coords, m = m)
  p <- ncol(Y)
  nr <- length(ranges)
  z_mat <- matrix(NA_real_, nrow = p, ncol = nr)
  u_mat <- matrix(NA_real_, nrow = p, ncol = nr)
  lambda_mat <- matrix(NA_real_, nrow = p, ncol = nr)
  colnames(z_mat) <- paste0("range_", seq_len(nr))
  colnames(u_mat) <- colnames(z_mat)
  colnames(lambda_mat) <- colnames(z_mat)
  rownames(z_mat) <- colnames(Y)
  rownames(u_mat) <- colnames(Y)
  rownames(lambda_mat) <- colnames(Y)

  for (a in seq_along(ranges)) {
    K <- spapath_sparse_kernel_from_knn(
      nn = nn,
      n = nrow(Y),
      range = ranges[a],
      kernel = kernel
    )
    moments <- .spapath_kernel_moments(K, X)
    KE <- as.matrix(K %*% E)
    quad <- colSums(E * KE)
    U <- 0.5 * (quad - moments$trace_mk)
    var_U <- 0.5 * moments$trace_mkmk
    z <- U / sqrt(pmax(var_U, .Machine$double.eps))
    lambda_hat <- pmax(0, (quad - moments$trace_mk) /
      pmax(moments$trace_mkmk, .Machine$double.eps))

    z_mat[, a] <- z
    u_mat[, a] <- U
    lambda_mat[, a] <- lambda_hat
  }

  best_index <- max.col(z_mat, ties.method = "first")
  z_best <- z_mat[cbind(seq_len(p), best_index)]
  u_best <- u_mat[cbind(seq_len(p), best_index)]
  lambda_best <- lambda_mat[cbind(seq_len(p), best_index)]
  spve_gene <- lambda_best / (1 + lambda_best)

  feature_scores <- data.frame(
    feature = colnames(Y),
    z = as.numeric(z_best),
    u = as.numeric(u_best),
    lambda = as.numeric(lambda_best),
    spve = as.numeric(spve_gene),
    best_range = ranges[best_index],
    stringsAsFactors = FALSE
  )

  list(
    feature_scores = feature_scores,
    z_by_range = z_mat,
    u_by_range = u_mat,
    lambda_by_range = lambda_mat,
    residuals = resid,
    sigma2 = sigma2,
    X = X,
    ranges = ranges
  )
}

spapath_pathway_tests <- function(scores,
                                  pathways,
                                  alpha_grid = c(0.5, 1, 2, 4, 8),
                                  topk_grid = NULL,
                                  n_sim = 2000,
                                  min_pathway_size = 2,
                                  max_pathway_size = 500) {
  feature_scores <- scores$feature_scores
  feature_index <- stats::setNames(seq_len(nrow(feature_scores)), feature_scores$feature)
  out <- vector("list", length(pathways))

  for (ii in seq_along(pathways)) {
    pname <- names(pathways)[ii]
    genes <- intersect(pathways[[ii]], feature_scores$feature)
    genes <- unique(genes)
    q <- length(genes)

    if (q < min_pathway_size || q > max_pathway_size) {
      out[[ii]] <- .spapath_empty_pathway_result(pname, q, "size_filter")
      next
    }

    idx <- feature_index[genes]
    resid_g <- scores$residuals[, idx, drop = FALSE]
    R <- .spapath_safe_cor(resid_g)
    q_eff <- .spapath_effective_size(R)
    Rz <- R^2
    diag(Rz) <- 1
    Rz <- .spapath_near_psd_cor(Rz)

    if (is.null(topk_grid)) {
      topk <- unique(pmin(q, c(1, 2, 5, 10, floor(sqrt(q)), q)))
    } else {
      topk <- unique(pmin(q, topk_grid))
    }
    topk <- topk[topk >= 1]

    z_null <- .spapath_sim_mvn(n_sim, Rz)
    null_stats <- .spapath_adaptive_stats_matrix(z_null, alpha_grid, topk)

    range_p <- numeric(ncol(scores$z_by_range))
    range_soft_p <- vector("list", ncol(scores$z_by_range))
    range_topk_p <- vector("list", ncol(scores$z_by_range))

    for (rr in seq_len(ncol(scores$z_by_range))) {
      z_rr <- scores$z_by_range[idx, rr]
      obs <- .spapath_adaptive_stats(z_rr, alpha_grid = alpha_grid, topk_grid = topk)
      p_soft <- .spapath_empirical_p(obs$soft, null_stats$soft)
      p_topk <- .spapath_empirical_p(obs$topk, null_stats$topk)
      range_soft_p[[rr]] <- p_soft
      range_topk_p[[rr]] <- p_topk
      range_p[rr] <- spapath_acat(c(p_soft, p_topk))
    }

    p_adapt <- spapath_acat(range_p)
    best_range_index <- which.min(range_p)
    best_soft <- which.min(range_soft_p[[best_range_index]])
    best_topk <- which.min(range_topk_p[[best_range_index]])
    alpha_hat <- alpha_grid[best_soft]
    z <- scores$z_by_range[idx, best_range_index]
    lambda <- scores$lambda_by_range[idx, best_range_index]
    spve <- lambda / (1 + lambda)
    driver <- .spapath_driver_genes(
      genes = genes,
      z = z,
      spve = spve,
      alpha = alpha_hat
    )

    e_spve_any <- mean(spve, na.rm = TRUE)
    e_spve_driver <- sum(driver$weights * driver$spve)
    e_spve_shrunk <- e_spve_any * q_eff / (q_eff + 5)

    out[[ii]] <- data.frame(
      pathway = pname,
      p_value = p_adapt,
      fdr = NA_real_,
      q = q,
      q_eff = q_eff,
      eSPVE_any = e_spve_any,
      eSPVE_any_shrunk = e_spve_shrunk,
      eSPVE_driver = e_spve_driver,
      best_alpha = alpha_hat,
      best_topk = topk[best_topk],
      best_range = scores$ranges[best_range_index],
      best_range_p = range_p[best_range_index],
      best_soft_p = min(range_soft_p[[best_range_index]]),
      best_topk_p = min(range_topk_p[[best_range_index]]),
      driver_genes = I(list(driver$genes)),
      driver_weights = I(list(driver$weights)),
      driver_scores = I(list(driver$z)),
      status = "tested",
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

spapath_default_ranges <- function(coords, multipliers = c(0.05, 0.15, 0.35)) {
  coords <- as.matrix(coords)
  bbox <- apply(coords, 2, range)
  diameter <- sqrt(sum((bbox[2, ] - bbox[1, ])^2))
  if (!is.finite(diameter) || diameter <= 0) {
    stop("coords must span a positive spatial domain")
  }
  diameter * multipliers
}

spapath_knn <- function(coords, m = 20) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  m <- min(as.integer(m), n - 1L)
  if (m < 1) {
    stop("m must be at least 1 and smaller than nrow(coords)")
  }

  if (requireNamespace("FNN", quietly = TRUE)) {
    fit <- FNN::get.knn(coords, k = m)
    return(list(index = fit$nn.index, dist = fit$nn.dist))
  }

  if (n > 5000) {
    stop("FNN is required for n > 5000. Install FNN or reduce n.")
  }

  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  index <- t(apply(D, 1, function(x) order(x)[seq_len(m)]))
  dist <- matrix(NA_real_, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    dist[i, ] <- D[i, index[i, ]]
  }
  list(index = index, dist = dist)
}

spapath_sparse_kernel_from_knn <- function(nn,
                                           n,
                                           range,
                                           kernel = c("matern32", "gaussian", "exponential")) {
  kernel <- match.arg(kernel)
  from <- rep(seq_len(n), ncol(nn$index))
  to <- as.vector(nn$index)
  d <- as.vector(nn$dist)
  keep <- is.finite(d) & from != to
  from <- from[keep]
  to <- to[keep]
  d <- d[keep]

  a <- pmin(from, to)
  b <- pmax(from, to)
  w <- .spapath_corr(d, range = range, kernel = kernel)
  edge_df <- data.frame(a = a, b = b, w = w)
  edge_df <- stats::aggregate(w ~ a + b, data = edge_df, FUN = max)

  Matrix::sparseMatrix(
    i = c(edge_df$a, edge_df$b, seq_len(n)),
    j = c(edge_df$b, edge_df$a, seq_len(n)),
    x = c(edge_df$w, edge_df$w, rep(1, n)),
    dims = c(n, n)
  )
}

spapath_acat <- function(p, weights = NULL) {
  p <- as.numeric(p)
  ok <- is.finite(p) & !is.na(p)
  p <- p[ok]
  if (length(p) == 0) {
    return(NA_real_)
  }
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  if (is.null(weights)) {
    weights <- rep(1 / length(p), length(p))
  } else {
    weights <- as.numeric(weights[ok])
    weights <- weights / sum(weights)
  }
  stat <- sum(weights * tan((0.5 - p) * pi))
  pmin(pmax(0.5 - atan(stat) / pi, 0), 1)
}

print.spapath_result <- function(x, n = 10, ...) {
  cat("SpaPath result\n")
  cat("  pathways tested:", sum(x$results$status == "tested"), "\n")
  cat("  kernel:", x$kernel, "\n")
  cat("  neighbor size:", x$m, "\n")
  cat("  ranges:", paste(signif(x$ranges, 4), collapse = ", "), "\n\n")
  cols <- c("pathway", "p_value", "fdr", "q", "q_eff", "eSPVE_any",
            "eSPVE_driver", "best_range", "best_alpha", "best_topk", "status")
  print(utils::head(x$results[, cols, drop = FALSE], n), row.names = FALSE)
  invisible(x)
}

.spapath_check_y <- function(Y) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (any(!is.finite(Y))) {
    stop("Y contains non-finite values")
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("feature_", seq_len(ncol(Y)))
  }
  Y
}

.spapath_check_coords <- function(coords, n) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (nrow(coords) != n) {
    stop("nrow(coords) must match nrow(Y)")
  }
  if (ncol(coords) < 2) {
    stop("coords must have at least two columns")
  }
  if (any(!is.finite(coords))) {
    stop("coords contains non-finite values")
  }
  coords
}

.spapath_prepare_design <- function(X, n) {
  if (is.null(X)) {
    X <- matrix(1, nrow = n, ncol = 1)
    colnames(X) <- "intercept"
    return(X)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != n) {
    stop("nrow(X) must match nrow(Y)")
  }
  if (any(!is.finite(X))) {
    stop("X contains non-finite values")
  }
  if (!any(apply(X, 2, function(x) all(abs(x - 1) < 1e-12)))) {
    X <- cbind(intercept = 1, X)
  }
  X
}

.spapath_check_pathways <- function(pathways, feature_names) {
  if (!is.list(pathways) || length(pathways) == 0) {
    stop("pathways must be a non-empty named list")
  }
  if (is.null(names(pathways)) || any(names(pathways) == "")) {
    names(pathways) <- paste0("pathway_", seq_along(pathways))
  }
  lapply(pathways, function(x) intersect(unique(as.character(x)), feature_names))
}

.spapath_residualize <- function(Y, X) {
  qrx <- qr(X)
  resid <- qr.resid(qrx, Y)
  rank_x <- qrx$rank
  df <- nrow(Y) - rank_x
  if (df <= 0) {
    stop("design matrix X has no residual degrees of freedom")
  }
  sigma2 <- colSums(resid^2) / df
  sigma2 <- pmax(sigma2, .Machine$double.eps)
  list(resid = resid, sigma2 = sigma2, df = df, rank = rank_x)
}

.spapath_kernel_moments <- function(K, X) {
  qrx <- qr(X)
  Q <- qr.Q(qrx)[, seq_len(qrx$rank), drop = FALSE]
  KQ <- as.matrix(K %*% Q)
  QKQ <- crossprod(Q, KQ)

  trace_k <- sum(Matrix::diag(K))
  trace_hk <- sum(diag(QKQ))
  trace_mk <- trace_k - trace_hk

  trace_k2 <- sum(K * K)
  trace_hk2 <- sum(KQ * KQ)
  trace_hkhk <- sum(QKQ * QKQ)
  trace_mkmk <- trace_k2 - 2 * trace_hk2 + trace_hkhk
  trace_mkmk <- max(as.numeric(trace_mkmk), .Machine$double.eps)

  list(trace_mk = as.numeric(trace_mk), trace_mkmk = trace_mkmk)
}

.spapath_corr <- function(d, range, kernel) {
  x <- pmax(d / range, 0)
  if (kernel == "matern32") {
    return((1 + sqrt(3) * x) * exp(-sqrt(3) * x))
  }
  if (kernel == "gaussian") {
    return(exp(-0.5 * x^2))
  }
  if (kernel == "exponential") {
    return(exp(-x))
  }
  stop("unknown kernel")
}

.spapath_safe_cor <- function(X) {
  X <- as.matrix(X)
  q <- ncol(X)
  if (q == 1) {
    return(matrix(1, 1, 1))
  }
  R <- suppressWarnings(stats::cor(X))
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

.spapath_effective_size <- function(R) {
  q <- ncol(R)
  denom <- sum(R^2)
  if (!is.finite(denom) || denom <= 0) {
    return(q)
  }
  q^2 / denom
}

.spapath_near_psd_cor <- function(R, eps = 1e-8) {
  R <- (R + t(R)) / 2
  ev <- eigen(R, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  R2 <- ev$vectors %*% (vals * t(ev$vectors))
  d <- sqrt(pmax(diag(R2), eps))
  R2 <- R2 / tcrossprod(d)
  R2 <- (R2 + t(R2)) / 2
  diag(R2) <- 1
  R2
}

.spapath_sim_mvn <- function(n_sim, Sigma) {
  q <- ncol(Sigma)
  if (q == 1) {
    return(matrix(stats::rnorm(n_sim), nrow = n_sim, ncol = 1))
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(MASS::mvrnorm(n = n_sim, mu = rep(0, q), Sigma = Sigma))
  }
  ev <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(ev$values, 0)
  Z <- matrix(stats::rnorm(n_sim * q), nrow = n_sim)
  Z %*% (ev$vectors %*% diag(sqrt(vals), q, q))
}

.spapath_adaptive_stats <- function(z, alpha_grid, topk_grid) {
  soft <- vapply(alpha_grid, function(a) .spapath_softmax(z, a), numeric(1))
  topk <- vapply(topk_grid, function(k) sum(sort(z, decreasing = TRUE)[seq_len(k)]),
                 numeric(1))
  names(soft) <- paste0("soft_", alpha_grid)
  names(topk) <- paste0("topk_", topk_grid)
  list(soft = soft, topk = topk)
}

.spapath_adaptive_stats_matrix <- function(Z, alpha_grid, topk_grid) {
  soft <- sapply(alpha_grid, function(a) {
    apply(Z, 1, .spapath_softmax, alpha = a)
  })
  topk <- sapply(topk_grid, function(k) {
    apply(Z, 1, function(z) sum(sort(z, decreasing = TRUE)[seq_len(k)]))
  })
  if (length(alpha_grid) == 1) {
    soft <- matrix(soft, ncol = 1)
  }
  if (length(topk_grid) == 1) {
    topk <- matrix(topk, ncol = 1)
  }
  colnames(soft) <- paste0("soft_", alpha_grid)
  colnames(topk) <- paste0("topk_", topk_grid)
  list(soft = soft, topk = topk)
}

.spapath_softmax <- function(z, alpha) {
  z <- as.numeric(z)
  zmax <- max(z)
  zmax + log(sum(exp(alpha * (z - zmax)))) / alpha
}

.spapath_empirical_p <- function(obs, null_mat) {
  obs <- as.numeric(obs)
  if (is.null(dim(null_mat))) {
    null_mat <- matrix(null_mat, ncol = 1)
  }
  p <- numeric(length(obs))
  for (j in seq_along(obs)) {
    p[j] <- (1 + sum(null_mat[, j] >= obs[j])) / (nrow(null_mat) + 1)
  }
  names(p) <- names(obs)
  p
}

.spapath_driver_genes <- function(genes, z, spve, alpha) {
  z_use <- z
  z_use[!is.finite(z_use)] <- min(z_use[is.finite(z_use)], 0)
  w <- exp(alpha * (z_use - max(z_use)))
  w[z <= 0] <- 0
  if (sum(w) <= 0) {
    w <- exp(alpha * (z_use - max(z_use)))
  }
  w <- w / sum(w)
  ord <- order(w, decreasing = TRUE)
  cumw <- cumsum(w[ord])
  keep <- ord[seq_len(max(1, which(cumw >= 0.8)[1]))]
  list(
    genes = genes[keep],
    weights = w[keep] / sum(w[keep]),
    z = z[keep],
    spve = spve[keep]
  )
}

.spapath_empty_pathway_result <- function(pathway, q, status) {
  data.frame(
    pathway = pathway,
    p_value = NA_real_,
    fdr = NA_real_,
    q = q,
    q_eff = NA_real_,
    eSPVE_any = NA_real_,
    eSPVE_any_shrunk = NA_real_,
    eSPVE_driver = NA_real_,
    best_alpha = NA_real_,
    best_topk = NA_integer_,
    best_range = NA_real_,
    best_range_p = NA_real_,
    best_soft_p = NA_real_,
    best_topk_p = NA_real_,
    driver_genes = I(list(character())),
    driver_weights = I(list(numeric())),
    driver_scores = I(list(numeric())),
    status = status,
    stringsAsFactors = FALSE
  )
}
