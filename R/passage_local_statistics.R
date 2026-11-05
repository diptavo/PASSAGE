# Pathway-local fitting is repeated for EVERY target/background set. The
# coordinate/design kernel cache is shared; no inner rotation tests are needed.
passage_make_local_gp_statistic <- function(Y, coords, X = NULL, method = "pca",
                                           K = NULL, variance_threshold = 0.95,
                                           kernel_cache = NULL,
                                           ranges = c(0.05, 0.1, 0.2, 0.4),
                                           kernel_method = "lowrank", n_features = 64L,
                                           aggregations = c("sum", "max"), top_k = 2L,
                                           fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
                                           seed = 1L) {
  Y <- passage_check_y(Y)
  if (is.null(rownames(Y))) rownames(Y) <- paste0("spot_", seq_len(nrow(Y)))
  coords <- passage_check_coords(coords, nrow(Y))
  if (is.null(rownames(coords))) rownames(coords) <- rownames(Y)
  if (!identical(rownames(coords), rownames(Y))) stop("coordinate spot IDs must exactly match Y")
  if (!is.null(X) && !is.null(rownames(X)) && !identical(rownames(X), rownames(Y)))
    stop("design spot IDs must exactly match Y")
  method <- match.arg(method, c("pca", "factor_analysis"))
  aggregations <- as.character(aggregations)
  if (!length(aggregations) || anyNA(aggregations) || anyDuplicated(aggregations) ||
      any(!aggregations %in% c("sum", "max", "topk"))) stop("aggregations must be unique sum/max/topk choices")
  if (length(top_k) != 1L || !is.finite(top_k) || top_k < 1 || top_k != as.integer(top_k))
    stop("top_k must be a positive integer")
  X <- passage_prepare_design(X, nrow(Y))
  rownames(X) <- rownames(Y)
  if (is.null(kernel_cache)) kernel_cache <- passage_prepare_gp_kernel(coords,
    X = X, ranges = ranges, kernel_method = kernel_method, n_features = n_features)
  else {
    if (!inherits(kernel_cache, "passage_gp_kernel") || !isTRUE(all.equal(coords, kernel_cache$coords)))
      stop("kernel_cache coordinates must exactly match coords")
    if (!isTRUE(all.equal(X, kernel_cache$X, check.attributes = TRUE)))
      stop("kernel_cache design must exactly match X")
    if (!missing(ranges) && !identical(sort(as.numeric(ranges)), kernel_cache$ranges))
      stop("ranges conflict with kernel_cache")
    if (!missing(kernel_method) && !identical(kernel_method, kernel_cache$kernel_method))
      stop("kernel_method conflicts with kernel_cache")
    if (!missing(n_features) && n_features != kernel_cache$n_features)
      stop("n_features conflicts with kernel_cache")
  }
  diag_store <- new.env(parent = emptyenv())
  diag_store$rows <- list()
  statistic <- function(idx) {
    if (is.numeric(idx) && (anyNA(idx) || any(!is.finite(idx)) ||
        any(idx != as.integer(idx)) || any(idx < 1L | idx > ncol(Y)))) stop("invalid pathway indices")
    idx <- passage_resolve_pathway(idx, colnames(Y))
    if (!length(idx)) stop("pathway has no genes present in Y")
    fit <- passage_fit_pathway(Y[, idx, drop = FALSE], seq_along(idx), X = X, method = method, K = K,
      variance_threshold = variance_threshold, seed = seed)
    detail <- passage_pathway_gp_statistic(fit, kernel_cache, fractions = fractions,
      aggregation = "sum", top_k = top_k, output = "details")
    per_pc <- detail$component_fits$statistic
    # Use the exact profile contributions computed once for all aggregations.
    if (is.null(per_pc)) stop("local GP profile did not provide component statistics")
    values <- vapply(aggregations, function(a) switch(a,
      sum = sum(per_pc), max = max(per_pc),
      topk = sum(sort(per_pc, decreasing = TRUE)[seq_len(min(top_k, length(per_pc)))])), numeric(1))
    names(values) <- paste0("local_gp_", aggregations,
      ifelse(aggregations == "topk", paste0(top_k), ""))
    row <- list(K = fit$K, explained_variance = fit$explained_variance)
    if (!length(diag_store$rows)) {
      row$component_fits <- detail$component_fits
      row$rank_rule <- fit$rank_rule
    }
    diag_store$rows[[length(diag_store$rows) + 1L]] <- row
    values
  }
  attr(statistic, "diagnostics") <- function() diag_store$rows
  attr(statistic, "reset_diagnostics") <- function() diag_store$rows <- list()
  attr(statistic, "configuration") <- list(estimation = "pathway", method = method,
    K = K, variance_threshold = variance_threshold, aggregations = aggregations,
    top_k = top_k, ranges = kernel_cache$ranges, fractions = fractions,
    kernel_method = kernel_cache$kernel_method, n_features = kernel_cache$n_features,
    design_rank = qr(X)$rank, calibration = "none; statistic only",
    limitation = "pathway-local competitive calibration remains research-stage")
  statistic
}
