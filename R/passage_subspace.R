# Scale-free trace and generalized largest root on pathway components.
passage_subspace_test <- function(Y, coords = NULL, X = NULL, spatial_kernel = NULL,
                                  kernel_range = 0.25, variance_threshold = 0.95,
                                  K = NULL, B = 9999L, seed = 1L, max_n = 3000L) {
  Y <- passage_check_y(Y)
  n <- nrow(Y)
  if (n > max_n) stop("dense reference implementation exceeds max_n; prespecify a spatial subsample")
  X <- passage_prepare_design(X, n)
  qx <- qr(X); df <- n - qx$rank
  if (qx$rank != ncol(X) || df < 3L) stop("invalid or saturated adjustment design")
  if (is.null(spatial_kernel)) {
    coords <- passage_check_coords(coords, n)
    if (length(kernel_range) != 1L || !is.finite(kernel_range) || kernel_range <= 0)
      stop("kernel_range must be positive")
    z <- sqrt(3) * as.matrix(stats::dist(coords)) / kernel_range
    spatial_kernel <- (1 + z) * exp(-z)
  }
  R <- as.matrix(spatial_kernel)
  if (!identical(dim(R), c(n, n)) || !all(is.finite(R)) || max(abs(R - t(R))) > 1e-8)
    stop("spatial_kernel must be a finite symmetric n by n matrix")
  R <- qr.resid(qx, t(qr.resid(qx, R)))
  R <- (R + t(R)) / 2
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(lambda) < -1e-7 * max(1, max(lambda))) stop("spatial_kernel must be positive semidefinite")
  lambda <- pmax(lambda[seq_len(df)], 0)
  if (sum(lambda) <= .Machine$double.eps) stop("kernel has no residual spatial variation")
  sy <- svd(qr.resid(qx, Y), nu = min(ncol(Y), df), nv = 0L)
  nr <- sum(sy$d > max(sy$d) * 1e-10)
  if (nr < 1L) stop("expression has no residual rank")
  if (is.null(K)) {
    if (!is.finite(variance_threshold) || variance_threshold <= 0 || variance_threshold > 1)
      stop("variance_threshold must be in (0,1]")
    K <- which(cumsum(sy$d^2) / sum(sy$d^2) >= variance_threshold)[1L]
    K <- min(K, nr)
  }
  if (length(K) != 1L || !is.finite(K) || K < 1 || K > nr || K != as.integer(K) || K >= df)
    stop("K must be a positive integer below residual df and at most expression rank")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("invalid B")
  U <- sy$u[, seq_len(K), drop = FALSE]
  stats_from <- function(A) c(trace = df * sum(diag(A)),
    largest_root = max(eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values))
  observed <- stats_from(crossprod(U, R %*% U))
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  null <- replicate(B, {
    Q <- qr.Q(qr(matrix(stats::rnorm(df * K), df, K)))
    stats_from(crossprod(Q, lambda * Q))
  })
  p <- (1 + rowSums(null >= observed)) / (B + 1)
  list(statistic = observed, p = p, p_omnibus_bonferroni = min(1, 2 * min(p)),
       null = null, K = K, residual_df = df, B = B, seed = seed,
       explained_variance = sum(sy$d[seq_len(K)]^2) / sum(sy$d^2),
       null_assumption = "iid Gaussian rows conditional on a fixed design; arbitrary gene covariance",
       role = "self-contained spatial subspace test, not a competitive pathway test")
}
