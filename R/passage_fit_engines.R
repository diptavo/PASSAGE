# Additional PASSAGE model-fitting engines for
#   Y_i = X_i B + A u(s_i) + e_i.
#
# These engines return the same "passage_engine" object shape as the original
# PCA/Vecchia engine, so existing score tests and effect-size metrics can reuse
# them without special-case code.

passage_factor_normalize <- function(V, min_sd = 1e-8) {
  V <- as.matrix(V)
  storage.mode(V) <- "double"
  if (ncol(V) == 0L) stop("V must contain at least one factor")
  V <- sweep(V, 2L, colMeans(V), "-")
  ss <- apply(V, 2L, stats::sd)
  keep <- is.finite(ss) & ss > min_sd
  if (!any(keep)) stop("all fitted factors are numerically degenerate")
  V <- V[, keep, drop = FALSE]
  ss <- ss[keep]
  V <- sweep(V, 2L, ss, "/")
  colnames(V) <- paste0("factor_", seq_len(ncol(V)))
  V
}

passage_refit_loadings <- function(R, V) {
  R <- as.matrix(R)
  V <- passage_factor_normalize(V)
  qrv <- qr(V)
  A <- t(qr.coef(qrv, R))
  A[!is.finite(A)] <- 0
  fitted_spatial <- V %*% t(A)
  D <- colMeans((R - fitted_spatial)^2)
  D <- pmax(D, .Machine$double.eps)
  rownames(A) <- colnames(R)
  colnames(A) <- colnames(V)
  colnames(fitted_spatial) <- colnames(R)
  list(A = A, V = V, fitted_spatial = fitted_spatial, D = D)
}

passage_engine_complete <- function(Y,
                                    coords,
                                    X,
                                    residuals,
                                    A,
                                    V,
                                    D,
                                    range_grid,
                                    m,
                                    kernel,
                                    ordering,
                                    engine_name,
                                    rank_info = NULL,
                                    fitted_spatial = NULL,
                                    extra = list(),
                                    verbose = TRUE) {
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  residuals <- as.matrix(residuals)
  A <- as.matrix(A)
  V <- as.matrix(V)
  K <- ncol(A)
  if (K != ncol(V)) stop("A and V must have the same number of factors")
  if (is.null(fitted_spatial)) fitted_spatial <- V %*% t(A)
  colnames(fitted_spatial) <- colnames(Y)
  rownames(A) <- colnames(Y)
  colnames(A) <- colnames(V) <- paste0("factor_", seq_len(K))

  if (is.null(range_grid)) range_grid <- passage_default_range_grid(coords)
  range_grid <- as.numeric(range_grid)
  if (any(!is.finite(range_grid)) || any(range_grid <= 0)) {
    stop("range_grid must contain positive finite values")
  }

  factor_fits <- vector("list", K)
  Q_list <- vector("list", K)
  K_score_list <- vector("list", K)
  vecchia_list <- vector("list", K)
  theta <- data.frame(
    sigma2 = rep(NA_real_, K),
    range = rep(NA_real_, K),
    effective_range = rep(NA_real_, K)
  )
  for (kk in seq_len(K)) {
    if (verbose) message(engine_name, ": fitting Vecchia GP parameters for factor ", kk, "/", K)
    ff <- passage_vecchia_fit_range_grid(
      v = V[, kk],
      coords = coords,
      range_grid = range_grid,
      m = m,
      kernel = kernel,
      ordering = ordering
    )
    factor_fits[[kk]] <- ff
    vecchia_list[[kk]] <- ff$vecchia
    Q_list[[kk]] <- ff$vecchia$Q
    K_score_list[[kk]] <- passage_sparse_covariance_from_vecchia(ff$vecchia)
    theta$sigma2[kk] <- ff$sigma2
    theta$range[kk] <- ff$range
    theta$effective_range[kk] <- passage_effective_range(ff$range, kernel)
  }
  rownames(theta) <- paste0("factor_", seq_len(K))

  out <- c(list(
    A = A,
    theta = theta,
    D = stats::setNames(pmax(as.numeric(D), .Machine$double.eps), colnames(Y)),
    V = V,
    Q = Q_list,
    K_score = K_score_list,
    vecchia = vecchia_list,
    factor_fits = factor_fits,
    X = X,
    residuals = residuals,
    fitted_spatial = fitted_spatial,
    gene_names = colnames(Y),
    coords = coords,
    center = TRUE,
    scale = FALSE,
    gene_center = rep(0, ncol(Y)),
    gene_scale = rep(1, ncol(Y)),
    kernel = kernel,
    m = m,
    range_grid = range_grid,
    rank_info = rank_info,
    N = nrow(Y),
    G = ncol(Y),
    K = K,
    engine = engine_name
  ), extra)
  class(out) <- c("passage_engine", "list")
  out
}

passage_svd_factor_fit <- function(R,
                                   F,
                                   K = 6L,
                                   rank_method = c("fixed", "variance"),
                                   variance_threshold = 0.90,
                                   min_K = 2L,
                                   max_K = NULL) {
  rank_method <- match.arg(rank_method)
  R <- as.matrix(R)
  F <- as.matrix(F)
  max_rank_possible <- min(ncol(R), nrow(R) - 1L, ncol(F), nrow(F))
  if (is.null(max_K)) max_K <- if (rank_method == "variance") min(50L, max_rank_possible) else as.integer(K)
  max_K <- min(as.integer(max_K), max_rank_possible)
  min_K <- min(as.integer(min_K), max_K)
  if (max_K < 1L) stop("K/max_K must be at least 1 after accounting for data dimensions")
  s <- passage_truncated_svd(F, max_K)
  rank_info <- passage_select_factor_rank(
    d = s$d,
    n = nrow(F),
    total_variance = sum(F^2) / max(1, nrow(F) - 1L),
    K = K,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    min_K = min_K,
    max_K = max_K
  )
  kk <- rank_info$K
  V0 <- s$u[, seq_len(kk), drop = FALSE] %*% diag(s$d[seq_len(kk)], kk, kk)
  fit <- passage_refit_loadings(R, V0)
  fit$rank_info <- rank_info
  fit
}

passage_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

passage_scale_columns <- function(Z, min_sd = 1e-8) {
  Z <- as.matrix(Z)
  Z <- sweep(Z, 2L, colMeans(Z), "-")
  ss <- apply(Z, 2L, stats::sd)
  keep <- is.finite(ss) & ss > min_sd
  if (!any(keep)) stop("spatial basis has no non-degenerate columns")
  Z <- Z[, keep, drop = FALSE]
  ss <- ss[keep]
  sweep(Z, 2L, ss, "/")
}

passage_spatial_basis_matrix <- function(coords,
                                         n_basis = 50L,
                                         basis = c("rbf", "polynomial", "rbf_poly"),
                                         knot_method = c("kmeans", "sample"),
                                         rbf_scale = NULL,
                                         include_polynomial = TRUE,
                                         seed = 20260526L) {
  basis <- match.arg(basis)
  knot_method <- match.arg(knot_method)
  coords <- passage_check_coords(coords)
  n <- nrow(coords)
  zc <- passage_scale_columns(coords)
  out <- list()

  if (basis %in% c("polynomial", "rbf_poly") || include_polynomial) {
    out$linear <- zc
    out$quadratic <- cbind(zc[, 1L]^2, zc[, 2L]^2, zc[, 1L] * zc[, 2L])
    colnames(out$quadratic) <- c("x2", "y2", "xy")
  }

  if (basis %in% c("rbf", "rbf_poly")) {
    nb <- min(as.integer(n_basis), n)
    centers <- passage_with_seed(seed, {
      if (knot_method == "kmeans" && nb < n) {
        tryCatch(stats::kmeans(zc, centers = nb, iter.max = 50)$centers,
                 error = function(e) zc[sample(seq_len(n), nb), , drop = FALSE])
      } else {
        zc[sample(seq_len(n), nb), , drop = FALSE]
      }
    })
    dc <- as.matrix(stats::dist(centers))
    if (is.null(rbf_scale)) {
      nn <- apply(dc + diag(Inf, nrow(dc)), 1L, min)
      rbf_scale <- stats::median(nn[is.finite(nn) & nn > 0])
      if (!is.finite(rbf_scale) || rbf_scale <= 0) {
        d_all <- as.matrix(stats::dist(zc))
        rbf_scale <- stats::median(d_all[d_all > 0])
      }
    }
    rbf_scale <- max(as.numeric(rbf_scale), sqrt(.Machine$double.eps))
    D2 <- matrix(0, nrow = n, ncol = nrow(centers))
    for (jj in seq_len(nrow(centers))) {
      D2[, jj] <- rowSums((zc - matrix(centers[jj, ], n, ncol(zc), byrow = TRUE))^2)
    }
    out$rbf <- exp(-0.5 * D2 / rbf_scale^2)
  }

  Z <- do.call(cbind, out)
  Z <- passage_scale_columns(Z)
  colnames(Z) <- paste0("basis_", seq_len(ncol(Z)))
  attr(Z, "basis") <- basis
  attr(Z, "n_basis_requested") <- n_basis
  attr(Z, "rbf_scale") <- rbf_scale
  Z
}

passage_residualize_basis <- function(Z, X) {
  X <- passage_prepare_design(X, nrow(Z), intercept = TRUE)
  Zr <- qr.resid(qr(X), as.matrix(Z))
  passage_scale_columns(Zr)
}

passage_ridge_fit_matrix <- function(Z, Y, lambda) {
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  q <- ncol(Z)
  lambda <- max(as.numeric(lambda), 0)
  solve(crossprod(Z) + diag(lambda + .Machine$double.eps, q), crossprod(Z, Y))
}

passage_choose_ridge_lambda <- function(Z, Y, lambda_grid = NULL) {
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  q <- ncol(Z)
  zz <- crossprod(Z)
  ev <- eigen((zz + t(zz)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  base <- mean(diag(zz))
  if (!is.finite(base) || base <= 0) base <- 1
  if (is.null(lambda_grid)) {
    lambda_grid <- base * c(0, exp(seq(log(1e-4), log(1e2), length.out = 12L)))
  }
  lambda_grid <- sort(unique(pmax(as.numeric(lambda_grid), 0)))
  rows <- vector("list", length(lambda_grid))
  fits <- vector("list", length(lambda_grid))
  for (ii in seq_along(lambda_grid)) {
    lam <- lambda_grid[[ii]]
    coef <- solve(zz + diag(lam + .Machine$double.eps, q), crossprod(Z, Y))
    F <- Z %*% coef
    rss <- sum((Y - F)^2)
    edf <- sum(ev / pmax(ev + lam + .Machine$double.eps, .Machine$double.eps))
    gcv <- rss / pmax((nrow(Z) - edf)^2, .Machine$double.eps)
    rows[[ii]] <- data.frame(lambda = lam, edf = edf, rss = rss, gcv = gcv)
    fits[[ii]] <- list(coef = coef, fitted = F)
  }
  grid <- do.call(rbind, rows)
  best <- which.min(grid$gcv)
  list(
    lambda = grid$lambda[[best]],
    coef = fits[[best]]$coef,
    fitted = fits[[best]]$fitted,
    grid = grid
  )
}

passage_fit_engine_spatial_basis <- function(Y,
                                             coords,
                                             X = NULL,
                                             K = 6L,
                                             rank_method = c("fixed", "variance"),
                                             variance_threshold = 0.90,
                                             max_K = NULL,
                                             min_K = 2L,
                                             n_basis = 50L,
                                             basis = c("rbf_poly", "rbf", "polynomial"),
                                             knot_method = c("kmeans", "sample"),
                                             lambda_grid = NULL,
                                             range_grid = NULL,
                                             m = 20L,
                                             kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                             ordering = c("coordinate", "none"),
                                             seed = 20260526L,
                                             verbose = TRUE) {
  kernel <- match.arg(kernel)
  ordering <- match.arg(ordering)
  basis <- match.arg(basis)
  knot_method <- match.arg(knot_method)
  rank_method <- match.arg(rank_method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  fit0 <- passage_residualize(Y, X)
  R <- fit0$resid
  if (verbose) message("PASSAGE spatial-basis engine: building residualized spatial basis")
  Z <- passage_spatial_basis_matrix(
    coords = coords,
    n_basis = n_basis,
    basis = basis,
    knot_method = knot_method,
    seed = seed
  )
  Zr <- passage_residualize_basis(Z, X)
  if (verbose) message("PASSAGE spatial-basis engine: selecting ridge penalty")
  ridge <- passage_choose_ridge_lambda(Zr, R, lambda_grid = lambda_grid)
  ff <- passage_svd_factor_fit(
    R = R,
    F = ridge$fitted,
    K = K,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    min_K = min_K,
    max_K = max_K
  )
  passage_engine_complete(
    Y = Y,
    coords = coords,
    X = X,
    residuals = R,
    A = ff$A,
    V = ff$V,
    D = ff$D,
    range_grid = range_grid,
    m = m,
    kernel = kernel,
    ordering = ordering,
    engine_name = "spatial_basis_ridge_v1",
    rank_info = ff$rank_info,
    fitted_spatial = ff$fitted_spatial,
    extra = list(
      basis_matrix = Zr,
      basis_spec = list(basis = basis, knot_method = knot_method, n_basis = n_basis),
      ridge_lambda = ridge$lambda,
      ridge_grid = ridge$grid,
      ridge_coef = ridge$coef
    ),
    verbose = verbose
  )
}

passage_fit_engine_smoothed_pca <- function(Y,
                                            coords,
                                            X = NULL,
                                            K = 6L,
                                            rank_method = c("fixed", "variance"),
                                            variance_threshold = 0.90,
                                            max_K = NULL,
                                            min_K = 2L,
                                            n_basis = 50L,
                                            basis = c("rbf_poly", "rbf", "polynomial"),
                                            knot_method = c("kmeans", "sample"),
                                            lambda_grid = NULL,
                                            range_grid = NULL,
                                            m = 20L,
                                            kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                            ordering = c("coordinate", "none"),
                                            seed = 20260526L,
                                            verbose = TRUE) {
  kernel <- match.arg(kernel)
  ordering <- match.arg(ordering)
  basis <- match.arg(basis)
  knot_method <- match.arg(knot_method)
  rank_method <- match.arg(rank_method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  fit0 <- passage_residualize(Y, X)
  R <- fit0$resid
  max_rank_possible <- min(ncol(Y), nrow(Y) - qr(X)$rank)
  if (is.null(max_K)) max_K <- if (rank_method == "variance") min(50L, max_rank_possible) else as.integer(K)
  max_K <- min(as.integer(max_K), max_rank_possible)
  if (verbose) message("PASSAGE smoothed-PCA engine: computing initial residual PCA")
  s <- passage_truncated_svd(R, max_K)
  rank_info <- passage_select_factor_rank(
    d = s$d,
    n = nrow(R),
    total_variance = sum(R^2) / max(1, nrow(R) - 1L),
    K = K,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    min_K = min_K,
    max_K = max_K
  )
  kk <- rank_info$K
  V0 <- s$u[, seq_len(kk), drop = FALSE] %*% diag(s$d[seq_len(kk)], kk, kk)
  if (verbose) message("PASSAGE smoothed-PCA engine: smoothing PCA scores over tissue coordinates")
  Z <- passage_spatial_basis_matrix(
    coords = coords,
    n_basis = n_basis,
    basis = basis,
    knot_method = knot_method,
    seed = seed
  )
  Zr <- passage_residualize_basis(Z, X)
  ridge <- passage_choose_ridge_lambda(Zr, V0, lambda_grid = lambda_grid)
  ff <- passage_refit_loadings(R, ridge$fitted)
  passage_engine_complete(
    Y = Y,
    coords = coords,
    X = X,
    residuals = R,
    A = ff$A,
    V = ff$V,
    D = ff$D,
    range_grid = range_grid,
    m = m,
    kernel = kernel,
    ordering = ordering,
    engine_name = "smoothed_pca_basis_v1",
    rank_info = rank_info,
    fitted_spatial = ff$fitted_spatial,
    extra = list(
      basis_matrix = Zr,
      basis_spec = list(basis = basis, knot_method = knot_method, n_basis = n_basis),
      ridge_lambda = ridge$lambda,
      ridge_grid = ridge$grid,
      ridge_coef = ridge$coef
    ),
    verbose = verbose
  )
}

passage_fit_engine_nmf <- function(Y,
                                   coords,
                                   X = NULL,
                                   K = 6L,
                                   n_iter = 300L,
                                   smooth_scores = TRUE,
                                   n_basis = 50L,
                                   basis = c("rbf_poly", "rbf", "polynomial"),
                                   knot_method = c("kmeans", "sample"),
                                   lambda_grid = NULL,
                                   range_grid = NULL,
                                   m = 20L,
                                   kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                   ordering = c("coordinate", "none"),
                                   seed = 20260526L,
                                   verbose = TRUE) {
  kernel <- match.arg(kernel)
  ordering <- match.arg(ordering)
  basis <- match.arg(basis)
  knot_method <- match.arg(knot_method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  fit0 <- passage_residualize(Y, X)
  R <- fit0$resid
  K <- min(as.integer(K), ncol(Y), nrow(Y) - qr(X)$rank)
  if (K < 1L) stop("K must be at least 1 after accounting for data dimensions")
  if (verbose) message("PASSAGE NMF engine: fitting shifted nonnegative factorization")
  Xnn <- R - min(R)
  Xnn <- pmax(Xnn, 0)
  nmf <- passage_nmf_mu(Xnn, K = K, n_iter = n_iter, seed = seed)
  V0 <- nmf$W
  ridge <- NULL
  if (smooth_scores) {
    if (verbose) message("PASSAGE NMF engine: smoothing NMF scores over tissue coordinates")
    Z <- passage_spatial_basis_matrix(
      coords = coords,
      n_basis = n_basis,
      basis = basis,
      knot_method = knot_method,
      seed = seed
    )
    Zr <- passage_residualize_basis(Z, X)
    ridge <- passage_choose_ridge_lambda(Zr, V0, lambda_grid = lambda_grid)
    V0 <- pmax(ridge$fitted, 0)
  }
  ff <- passage_refit_loadings(R, V0)
  rank_info <- list(
    K = ncol(ff$A),
    requested_K = K,
    rank_method = "fixed",
    cumulative_variance = NA_real_,
    threshold_reached = NA
  )
  extra <- list(nmf = nmf, smooth_scores = smooth_scores, n_iter = n_iter)
  if (!is.null(ridge)) {
    extra$ridge_lambda <- ridge$lambda
    extra$ridge_grid <- ridge$grid
    extra$ridge_coef <- ridge$coef
  }
  passage_engine_complete(
    Y = Y,
    coords = coords,
    X = X,
    residuals = R,
    A = ff$A,
    V = ff$V,
    D = ff$D,
    range_grid = range_grid,
    m = m,
    kernel = kernel,
    ordering = ordering,
    engine_name = "nmf_spatial_basis_v1",
    rank_info = rank_info,
    fitted_spatial = ff$fitted_spatial,
    extra = extra,
    verbose = verbose
  )
}

passage_initial_factor_fit <- function(R,
                                       K,
                                       rank_method = c("fixed", "variance"),
                                       variance_threshold = 0.90,
                                       min_K = 2L,
                                       max_K = NULL) {
  rank_method <- match.arg(rank_method)
  max_rank_possible <- min(ncol(R), nrow(R) - 1L)
  if (is.null(max_K)) max_K <- if (rank_method == "variance") min(50L, max_rank_possible) else as.integer(K)
  max_K <- min(as.integer(max_K), max_rank_possible)
  s <- passage_truncated_svd(R, max_K)
  rank_info <- passage_select_factor_rank(
    d = s$d,
    n = nrow(R),
    total_variance = sum(R^2) / max(1, nrow(R) - 1L),
    K = K,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    min_K = min_K,
    max_K = max_K
  )
  kk <- rank_info$K
  V0 <- s$u[, seq_len(kk), drop = FALSE] %*% diag(s$d[seq_len(kk)], kk, kk)
  ff <- passage_refit_loadings(R, V0)
  ff$rank_info <- rank_info
  ff
}

passage_fit_engine_alternating_gp <- function(Y,
                                              coords,
                                              X = NULL,
                                              K = 6L,
                                              rank_method = c("fixed", "variance"),
                                              variance_threshold = 0.90,
                                              max_K = NULL,
                                              min_K = 2L,
                                              n_iter = 4L,
                                              smooth_penalty = 1,
                                              range_grid = NULL,
                                              m = 20L,
                                              kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                              ordering = c("coordinate", "none"),
                                              refit_range_every = 1L,
                                              verbose = TRUE) {
  kernel <- match.arg(kernel)
  ordering <- match.arg(ordering)
  rank_method <- match.arg(rank_method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  fit0 <- passage_residualize(Y, X)
  R <- fit0$resid
  if (is.null(range_grid)) range_grid <- passage_default_range_grid(coords)
  if (verbose) message("PASSAGE alternating-GP engine: initializing factors from residual PCA")
  ff <- passage_initial_factor_fit(
    R = R,
    K = K,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    min_K = min_K,
    max_K = max_K
  )
  A <- ff$A
  V <- ff$V
  D <- ff$D
  K_fit <- ncol(A)
  I_n <- Matrix::Diagonal(nrow(Y))
  factor_fits <- vector("list", K_fit)
  refit_range_every <- max(1L, as.integer(refit_range_every))

  for (iter in seq_len(as.integer(n_iter))) {
    if (verbose) message("PASSAGE alternating-GP engine: iteration ", iter, "/", n_iter)
    if (iter == 1L || (iter - 1L) %% refit_range_every == 0L) {
      for (kk in seq_len(K_fit)) {
        factor_fits[[kk]] <- passage_vecchia_fit_range_grid(
          v = V[, kk],
          coords = coords,
          range_grid = range_grid,
          m = m,
          kernel = kernel,
          ordering = ordering
        )
      }
    }
    Dinv <- 1 / pmax(D, .Machine$double.eps)
    for (kk in seq_len(K_fit)) {
      if (K_fit == 1L) {
        R_without <- R
      } else {
        other <- setdiff(seq_len(K_fit), kk)
        R_without <- R - V[, other, drop = FALSE] %*% t(A[, other, drop = FALSE])
      }
      rhs <- as.numeric(R_without %*% (A[, kk] * Dinv))
      alpha <- sum(A[, kk]^2 * Dinv)
      if (!is.finite(alpha) || alpha <= .Machine$double.eps) next
      Q <- factor_fits[[kk]]$vecchia$Q
      lhs <- Matrix::forceSymmetric(alpha * I_n + smooth_penalty * Q + 1e-8 * I_n, uplo = "U")
      V[, kk] <- as.numeric(Matrix::solve(lhs, rhs))
    }
    ff_iter <- passage_refit_loadings(R, V)
    A <- ff_iter$A
    V <- ff_iter$V
    D <- ff_iter$D
    K_fit <- ncol(A)
  }
  passage_engine_complete(
    Y = Y,
    coords = coords,
    X = X,
    residuals = R,
    A = A,
    V = V,
    D = D,
    range_grid = range_grid,
    m = m,
    kernel = kernel,
    ordering = ordering,
    engine_name = "alternating_gp_vecchia_v1",
    rank_info = ff$rank_info,
    fitted_spatial = V %*% t(A),
    extra = list(n_iter = n_iter, smooth_penalty = smooth_penalty),
    verbose = verbose
  )
}

passage_fit_factor_engine <- function(Y,
                                      coords,
                                      X = NULL,
                                      method = c("pca", "spatial_basis", "smoothed_pca", "nmf", "alternating_gp"),
                                      ...) {
  method <- match.arg(method)
  switch(method,
    pca = passage_fit_engine_pca(Y = Y, coords = coords, X = X, ...),
    spatial_basis = passage_fit_engine_spatial_basis(Y = Y, coords = coords, X = X, ...),
    smoothed_pca = passage_fit_engine_smoothed_pca(Y = Y, coords = coords, X = X, ...),
    nmf = passage_fit_engine_nmf(Y = Y, coords = coords, X = X, ...),
    alternating_gp = passage_fit_engine_alternating_gp(Y = Y, coords = coords, X = X, ...)
  )
}

passage_compare_factor_engines <- function(Y,
                                           coords,
                                           X = NULL,
                                           methods = c("pca", "spatial_basis", "smoothed_pca", "nmf", "alternating_gp"),
                                           fit_args = list(),
                                           verbose = TRUE) {
  fits <- vector("list", length(methods))
  rows <- vector("list", length(methods))
  names(fits) <- methods
  for (ii in seq_along(methods)) {
    method <- methods[[ii]]
    if (verbose) message("PASSAGE engine comparison: fitting ", method)
    args <- c(list(Y = Y, coords = coords, X = X, method = method), fit_args)
    fit <- tryCatch(do.call(passage_fit_factor_engine, args), error = function(e) e)
    fits[[ii]] <- fit
    rows[[ii]] <- if (inherits(fit, "error")) {
      data.frame(method = method, engine = NA_character_, status = "error",
                 K = NA_integer_, rss = NA_real_, mean_residual_var = NA_real_,
                 mean_spatial_var = NA_real_, error_message = conditionMessage(fit))
    } else {
      spatial_var <- colMeans(fit$fitted_spatial^2)
      data.frame(method = method, engine = fit$engine, status = "fit",
                 K = fit$K, rss = sum((fit$residuals - fit$fitted_spatial)^2),
                 mean_residual_var = mean(fit$D),
                 mean_spatial_var = mean(spatial_var),
                 error_message = NA_character_)
    }
  }
  list(summary = do.call(rbind, rows), fits = fits)
}
