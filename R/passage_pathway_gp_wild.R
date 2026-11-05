# Heteroskedasticity-aware residual bootstrap. One multiplier per spot acts on
# the whole gene vector, retaining within-spot gene covariance. Unlike Gaussian
# rotations, this is an approximate residual-bootstrap procedure, not exact.
passage_test_pathway_gp_wild <- function(Y, pathway, coords = NULL, X = NULL,
                                         method = c("pca", "factor_analysis"),
                                         K = NULL, variance_threshold = 0.95,
                                         kernel_cache = NULL,
                                         ranges = c(0.05, 0.1, 0.2, 0.4),
                                         fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
                                         kernel_method = c("lowrank", "dense"),
                                         n_features = 64L, B = 999L, seed = 1L,
                                         aggregation = c("sum", "max", "topk"),
                                         top_k = 2L, leverage = c("hc2", "none"),
                                         keep_null = TRUE) {
  started <- proc.time()[[3L]]
  supplied_kernel_method <- !missing(kernel_method)
  method <- match.arg(method)
  kernel_method <- match.arg(kernel_method)
  aggregation <- match.arg(aggregation)
  leverage <- match.arg(leverage)
  B <- .passage_gp_integer(B, "B")
  seed <- .passage_gp_integer(seed, "seed", 0L)
  top_k <- .passage_gp_integer(top_k, "top_k")
  if (!is.logical(keep_null) || length(keep_null) != 1L || is.na(keep_null)) stop("keep_null must be TRUE or FALSE")
  Y <- passage_check_y(Y)
  fit <- passage_fit_pathway(Y, pathway, X = X, method = method,
    K = K, variance_threshold = variance_threshold, seed = seed)
  Y <- Y[, fit$gene_names, drop = FALSE]
  rownames(Y) <- rownames(fit$scores)
  X <- fit$X
  if (is.null(kernel_cache)) {
    if (is.null(coords)) stop("supply coords or kernel_cache")
    if (is.null(rownames(coords))) rownames(coords) <- rownames(Y)
    kernel_cache <- passage_prepare_gp_kernel(coords, X = X, ranges = ranges,
      kernel_method = kernel_method, n_features = n_features)
  } else {
    if (!is.null(coords) && !isTRUE(all.equal(as.matrix(coords), kernel_cache$coords)))
      stop("coords conflict with kernel_cache")
    if (!missing(ranges) && !identical(sort(as.numeric(ranges)), kernel_cache$ranges))
      stop("ranges conflict with kernel_cache")
    if (supplied_kernel_method && !identical(kernel_method, kernel_cache$kernel_method))
      stop("kernel_method conflicts with kernel_cache")
    if (!missing(n_features) && n_features != kernel_cache$n_features)
      stop("n_features conflicts with kernel_cache")
  }
  observed <- passage_pathway_gp_statistic(fit, kernel_cache, fractions,
    aggregation, top_k, output = "details")
  qx <- qr(X)
  residuals <- qr.resid(qx, Y)
  fitted <- Y - residuals
  h <- rowSums(qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]^2)
  if (any(h >= 1 - 1e-8)) stop("near-unit leverage: residual wild bootstrap is unsupported")
  if (leverage == "hc2") residuals <- residuals / sqrt(1 - h)
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  null <- numeric(B)
  ranks <- integer(B)
  for (b in seq_len(B)) {
    multiplier <- 2L * sample.int(2L, nrow(Y), replace = TRUE) - 3L
    simulated <- fitted + residuals * multiplier
    # Repeat adjustment, PCA/FA, rank selection, and all GP parameter searches.
    fit_b <- passage_fit_pathway(simulated, colnames(Y), X = X, method = method,
      K = K, variance_threshold = variance_threshold, seed = seed)
    null[b] <- passage_pathway_gp_statistic(fit_b, kernel_cache, fractions,
      aggregation, top_k)
    ranks[b] <- fit_b$K
  }
  exceedances <- sum(null >= observed$statistic - 1e-10 * max(1, abs(observed$statistic)))
  list(statistic = observed$statistic, p = (1 + exceedances) / (B + 1),
    B = B, seed = seed, exceedances = exceedances, mc_floor = 1 / (B + 1),
    component_fits = observed$component_fits, estimation = observed$estimation,
    aggregation = aggregation, top_k = observed$top_k, top_k_used = observed$top_k_used,
    null = if (keep_null) null else NULL, null_ranks = if (keep_null) ranks else NULL,
    method = "gp_profile_wild_refit", leverage = leverage, max_leverage = max(h),
    calibration = "approximate Rademacher residual wild bootstrap; entire pathway estimation and GP search refitted",
    null_assumption = paste("Independent, mean-zero residual spot vectors after the specified mean design;",
      "spot variances and within-spot gene covariance may differ. Finite-sample validity is not guaranteed.",
      "Residual spatial dependence, misspecified means, high leverage and estimated covariates require validation."),
    role = "self-contained spatial covariance test; not competitive",
    elapsed_seconds = proc.time()[[3L]] - started)
}
