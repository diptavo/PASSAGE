# PASSAGE Layer 1: fast PCA two-stage spatial factor engine.

passage_fit_engine_pca <- function(Y,
                                   coords,
                                   X = NULL,
                                   K = 6L,
                                   rank_method = c("fixed", "variance"),
                                   variance_threshold = 0.90,
                                   max_K = NULL,
                                   min_K = 2L,
                                   factor_method = c("pca", "sparse_pca", "varimax", "ica", "nmf"),
                                   sparse_top_frac = 0.05,
                                   sparse_min_loadings = 10L,
                                   m = 20L,
                                   range_grid = NULL,
                                   kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                   ordering = c("coordinate", "none"),
                                   center = TRUE,
                                   scale = FALSE,
                                   verbose = TRUE) {
  kernel <- match.arg(kernel)
  rank_method <- match.arg(rank_method)
  factor_method <- match.arg(factor_method)
  ordering <- match.arg(ordering)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  max_rank_possible <- min(ncol(Y), nrow(Y) - qr(X)$rank)
  if (is.null(max_K)) max_K <- if (rank_method == "variance") min(50L, max_rank_possible) else as.integer(K)
  max_K <- min(as.integer(max_K), max_rank_possible)
  min_K <- min(as.integer(min_K), max_K)
  if (max_K < 1L) stop("K/max_K must be at least 1 after accounting for data dimensions")
  if (is.null(range_grid)) {
    range_grid <- passage_default_range_grid(coords)
  }
  range_grid <- as.numeric(range_grid)
  if (any(!is.finite(range_grid)) || any(range_grid <= 0)) {
    stop("range_grid must contain positive finite values")
  }

  if (verbose) {
    message("PASSAGE PCA engine: residualizing expression")
  }
  fit0 <- passage_residualize(Y, X)
  R <- fit0$resid
  gene_center <- rep(0, ncol(R))
  gene_scale <- rep(1, ncol(R))
  if (center) {
    gene_center <- colMeans(R)
    R <- sweep(R, 2L, gene_center, "-")
  }
  if (scale) {
    gene_scale <- apply(R, 2L, stats::sd)
    gene_scale <- pmax(gene_scale, sqrt(.Machine$double.eps))
    R <- sweep(R, 2L, gene_scale, "/")
  }

  if (verbose) {
    message("PASSAGE PCA engine: computing low-rank loading basis")
  }
  svd_rank <- if (rank_method == "variance") max_K else min(as.integer(K), max_K)
  s <- passage_truncated_svd(R, svd_rank)
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
  K <- rank_info$K
  if (verbose) {
    message(
      "PASSAGE PCA engine: selected K=", K,
      " (method=", rank_method,
      ", cumulative variance=", sprintf("%.3f", rank_info$cumulative_variance),
      if (!rank_info$threshold_reached) ", threshold not reached before max_K" else "",
      "); factor_method=", factor_method
    )
  }
  eig <- (s$d[seq_len(K)]^2) / max(1, nrow(R) - 1L)
  basis <- passage_low_rank_basis(
    R, s, K, eig, factor_method,
    sparse_top_frac = sparse_top_frac,
    sparse_min_loadings = sparse_min_loadings
  )
  A_hat <- basis$A
  V_scores <- basis$V
  fitted_spatial <- V_scores %*% t(A_hat)
  D_hat <- colMeans((R - fitted_spatial)^2)
  D_hat <- pmax(D_hat, .Machine$double.eps)

  if (verbose) {
    message("PASSAGE PCA engine: fitting Vecchia GP parameters per factor")
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
  for (k in seq_len(K)) {
    ff <- passage_vecchia_fit_range_grid(
      v = V_scores[, k],
      coords = coords,
      range_grid = range_grid,
      m = m,
      kernel = kernel,
      ordering = ordering
    )
    factor_fits[[k]] <- ff
    vecchia_list[[k]] <- ff$vecchia
    Q_list[[k]] <- ff$vecchia$Q
    K_score_list[[k]] <- passage_sparse_covariance_from_vecchia(ff$vecchia)
    theta$sigma2[k] <- ff$sigma2
    theta$range[k] <- ff$range
    theta$effective_range[k] <- passage_effective_range(ff$range, kernel)
  }

  out <- list(
    A = A_hat,
    theta = theta,
    D = D_hat,
    V = V_scores,
    Q = Q_list,
    K_score = K_score_list,
    vecchia = vecchia_list,
    factor_fits = factor_fits,
    X = X,
    residuals = R,
    gene_names = colnames(Y),
    coords = coords,
    center = center,
    scale = scale,
    gene_center = gene_center,
    gene_scale = gene_scale,
    kernel = kernel,
    m = m,
    range_grid = range_grid,
    rank_method = rank_method,
    variance_threshold = variance_threshold,
    rank_info = rank_info,
    factor_method = factor_method,
    sparse_top_frac = sparse_top_frac,
    sparse_min_loadings = sparse_min_loadings,
    N = nrow(Y),
    G = ncol(Y),
    K = K,
    engine = "pca_vecchia_v1"
  )
  class(out) <- c("passage_engine", "list")
  out
}

print.passage_engine <- function(x, ...) {
  cat("PASSAGE engine\n")
  cat("  type:", x$engine, "\n")
  cat("  N:", x$N, " G:", x$G, " K:", x$K, "\n")
  cat("  kernel:", x$kernel, " m:", x$m, "\n")
  print(x$theta)
  invisible(x)
}

passage_truncated_svd <- function(X, K) {
  if (requireNamespace("irlba", quietly = TRUE)) {
    fit <- irlba::irlba(X, nv = K, nu = K)
    return(list(d = fit$d, u = fit$u, v = fit$v))
  }
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    fit <- RSpectra::svds(X, k = K, nu = K, nv = K)
    return(list(d = fit$d, u = fit$u, v = fit$v))
  }
  svd(X, nu = K, nv = K)
}

passage_select_factor_rank <- function(d,
                                       n,
                                       total_variance = NULL,
                                       K,
                                       rank_method = c("fixed", "variance"),
                                       variance_threshold = 0.90,
                                       min_K = 2L,
                                       max_K = length(d)) {
  rank_method <- match.arg(rank_method)
  d <- as.numeric(d)
  eig <- d^2 / max(1, n - 1L)
  if (is.null(total_variance) || !is.finite(total_variance) || total_variance <= 0) {
    total_variance <- sum(eig)
  }
  cum <- cumsum(eig) / max(total_variance, .Machine$double.eps)
  if (rank_method == "fixed") {
    kk <- min(as.integer(K), max_K, length(d))
  } else {
    hit <- which(cum >= variance_threshold)
    kk <- if (length(hit)) hit[[1L]] else min(max_K, length(d))
    kk <- min(kk, max_K, length(d))
    kk <- max(kk, min_K)
  }
  list(
    K = kk,
    requested_K = K,
    max_K = max_K,
    min_K = min_K,
    variance_threshold = variance_threshold,
    cumulative_variance = cum[[kk]],
    threshold_reached = cum[[kk]] >= variance_threshold,
    singular_values = d,
    cumulative_variance_by_rank = cum
  )
}

passage_low_rank_basis <- function(R,
                                   s,
                                   K,
                                   eig,
                                   factor_method = c("pca", "sparse_pca", "varimax", "ica", "nmf"),
                                   sparse_top_frac = 0.05,
                                   sparse_min_loadings = 10L) {
  factor_method <- match.arg(factor_method)
  A_pca <- s$v[, seq_len(K), drop = FALSE] %*% diag(sqrt(pmax(eig, 0)), K, K)
  V_pca <- s$u[, seq_len(K), drop = FALSE] %*% diag(s$d[seq_len(K)], K, K)
  if (factor_method == "pca") {
    return(list(A = A_pca, V = V_pca, method = "pca"))
  }
  if (factor_method == "sparse_pca") {
    A <- A_pca
    n_keep <- max(as.integer(sparse_min_loadings), ceiling(nrow(A) * sparse_top_frac))
    n_keep <- min(nrow(A), n_keep)
    for (k in seq_len(ncol(A))) {
      keep <- order(abs(A[, k]), decreasing = TRUE)[seq_len(n_keep)]
      drop <- setdiff(seq_len(nrow(A)), keep)
      A[drop, k] <- 0
    }
    gram <- crossprod(A)
    gram <- gram + diag(.Machine$double.eps * max(diag(gram), 1), nrow(gram))
    V <- R %*% A %*% solve(gram)
    V <- scale(V)
    A <- t(qr.coef(qr(V), R))
    thresh <- apply(abs(A), 2L, function(x) sort(x, decreasing = TRUE)[min(n_keep, length(x))])
    A[sweep(abs(A), 2L, thresh, "<")] <- 0
    return(list(A = as.matrix(A), V = as.matrix(V), method = "sparse_pca_threshold"))
  }
  if (factor_method == "varimax") {
    rot <- stats::varimax(A_pca)
    Tmat <- as.matrix(rot$rotmat)
    A <- as.matrix(rot$loadings)
    V <- V_pca %*% Tmat
    return(list(A = A, V = V, method = "varimax"))
  }
  if (factor_method == "ica") {
    if (requireNamespace("fastICA", quietly = TRUE)) {
      fit <- fastICA::fastICA(V_pca, n.comp = K, method = "C", verbose = FALSE)
      V <- scale(fit$S)
      A <- t(qr.coef(qr(V), R))
      return(list(A = A, V = as.matrix(V), method = "ica"))
    }
    warning("fastICA is not installed; using varimax-rotated PCA as ICA fallback")
    return(passage_low_rank_basis(R, s, K, eig, factor_method = "varimax"))
  }
  if (factor_method == "nmf") {
    X <- R - min(R)
    X <- pmax(X, 0)
    nmf <- passage_nmf_mu(X, K = K, n_iter = 200L)
    V <- scale(nmf$W)
    A <- t(qr.coef(qr(V), R))
    return(list(A = A, V = as.matrix(V), method = "nmf_mu_shifted"))
  }
  stop("unknown factor_method")
}

passage_nmf_mu <- function(X, K, n_iter = 200L, seed = 20260523L, eps = 1e-8) {
  set.seed(seed)
  X <- as.matrix(X)
  W <- matrix(stats::runif(nrow(X) * K, min = 0.01, max = 1), nrow(X), K)
  H <- matrix(stats::runif(K * ncol(X), min = 0.01, max = 1), K, ncol(X))
  for (ii in seq_len(n_iter)) {
    H <- H * (crossprod(W, X) / pmax(crossprod(W, W) %*% H, eps))
    W <- W * ((X %*% t(H)) / pmax(W %*% (H %*% t(H)), eps))
  }
  list(W = W, H = H)
}
