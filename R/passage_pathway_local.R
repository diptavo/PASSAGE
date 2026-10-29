# Pathway-local estimation is independent of coordinates and the global engine.
passage_fit_pathway <- function(Y, pathway, X = NULL,
                                method = c("pca", "factor_analysis"),
                                K = NULL, variance_threshold = 0.95,
                                seed = 1L) {
  method <- match.arg(method)
  Y <- passage_check_y(Y)
  if (nrow(Y) < 5L || ncol(Y) < 2L || anyDuplicated(colnames(Y)) ||
      anyNA(colnames(Y)) || any(!nzchar(colnames(Y))))
    stop("Y needs at least five spots and unique nonempty gene names")
  if (is.null(rownames(Y))) rownames(Y) <- paste0("spot_", seq_len(nrow(Y)))
  if (anyDuplicated(rownames(Y)) || anyNA(rownames(Y)) || any(!nzchar(rownames(Y))))
    stop("Y needs unique nonempty spot names")
  if (is.numeric(pathway) && (anyNA(pathway) || any(!is.finite(pathway)) ||
      any(pathway != as.integer(pathway)) || any(pathway < 1L | pathway > ncol(Y))))
    stop("pathway indices must be integers within Y")
  if (!is.numeric(pathway) && (!is.character(pathway) || anyNA(pathway) || any(!nzchar(pathway))))
    stop("pathway must contain gene names or indices")
  P <- passage_resolve_pathway(pathway, colnames(Y))
  if (length(P) < 2L) stop("fewer than two pathway genes are present")
  missing <- if (is.character(pathway)) setdiff(pathway, colnames(Y)) else character()
  if (!is.null(X) && !is.null(rownames(X)) && !identical(rownames(X), rownames(Y)))
    stop("X spot IDs must exactly match Y in order")
  X <- passage_prepare_design(X, nrow(Y))
  rownames(X) <- rownames(Y)
  qx <- qr(X)
  df <- nrow(Y) - qx$rank
  if (qx$rank != ncol(X) || df < 3L) stop("X must be full rank with at least three residual degrees of freedom")
  residuals <- qr.resid(qx, Y[, P, drop = FALSE])
  variances <- colSums(residuals^2) / df
  keep <- variances > max(variances) * 1e-12 & variances > 0
  dropped <- colnames(residuals)[!keep]
  residuals <- residuals[, keep, drop = FALSE]
  if (ncol(residuals) < 2L) stop("fewer than two nonconstant pathway genes after adjustment")
  if (!is.null(K) && (length(K) != 1L || !is.finite(K) || K < 1 || K != as.integer(K)))
    stop("K must be a positive integer")
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  diagnostics <- list()
  if (method == "pca") {
    s <- svd(residuals, nu = min(ncol(residuals), df), nv = min(ncol(residuals), df))
    nr <- sum(s$d > max(s$d) * 1e-10)
    if (is.null(K)) {
      if (length(variance_threshold) != 1L || !is.finite(variance_threshold) ||
          variance_threshold <= 0 || variance_threshold > 1)
        stop("variance_threshold must be in (0,1]")
      K <- min(which(cumsum(s$d^2) / sum(s$d^2) >= variance_threshold)[1L], nr)
      rank_rule <- "pathway PCA cumulative variance"
    } else rank_rule <- "prespecified rank"
    if (K > nr || K >= df) stop("K must not exceed expression rank and must be below residual df")
    scores <- sweep(s$u[, seq_len(K), drop = FALSE], 2L, s$d[seq_len(K)], "*")
    loadings <- s$v[, seq_len(K), drop = FALSE]
    explained_variance <- sum(s$d[seq_len(K)]^2) / sum(s$d^2)
    diagnostics$singular_values <- s$d
    diagnostics$gene_scaling <- "none; covariance PCA without gene variance scaling"
  } else {
    if (is.null(K)) stop("factor_analysis requires a prespecified K; PCA variance retention is not an FA rank rule")
    if (K >= df) stop("K must be below residual df")
    warnings <- character()
    fa <- withCallingHandlers(
      stats::factanal(residuals, factors = K, scores = "regression", rotation = "none"),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    # factanal uses the gene correlation matrix; convert loadings back to expression units.
    gene_sd <- apply(residuals, 2L, stats::sd)
    scores <- as.matrix(fa$scores)
    # Drop factanal's special 'loadings' class so exported tables remain rectangular.
    loadings <- sweep(unclass(fa$loadings), 1L, gene_sd, "*")
    if (qr(scores)$rank != K) stop("factor-analysis scores are rank deficient")
    rank_rule <- "prespecified rank"
    explained_variance <- NA_real_
    diagnostics <- list(gene_scaling = "gene correlation matrix; regression factor scores",
      uniqueness = fa$uniquenesses, noise_variance = fa$uniquenesses * gene_sd^2,
      boundary_uniqueness = any(fa$uniquenesses <= 0.005 * (1 + 1e-6)), warnings = warnings,
      note = "FA extension to the PCA baseline; factanal model-fit p-values are not spatial tests")
  }
  rownames(scores) <- rownames(Y)
  rownames(loadings) <- colnames(residuals)
  colnames(scores) <- colnames(loadings) <- paste0("component_", seq_len(K))
  reconstruction <- scores %*% t(loadings)
  noise <- colSums((residuals - reconstruction)^2) / df
  out <- list(scores = scores, loadings = loadings, X = X, residual_df = df,
    gene_names = colnames(residuals), missing_genes = missing, dropped_genes = dropped,
    K = as.integer(K), method = method, rank_rule = rank_rule,
    explained_variance = explained_variance, residual_variance = noise,
    reconstruction_fraction = 1 - sum((residuals - reconstruction)^2) / sum(residuals^2),
    scores_source = "unsmoothed_expression", estimation_scope = "pathway_only",
    seed = seed, diagnostics = diagnostics)
  class(out) <- c("passage_pathway_fit", "list")
  out
}

# Test the estimated pathway subspace. GP fits are deliberately a separate object.
passage_test_pathway_fit <- function(fit, coords = NULL, spatial_kernel = NULL,
                                     kernel_range = 0.25, B = 9999L,
                                     seed = 1L, max_n = 3000L) {
  if (!inherits(fit, "passage_pathway_fit") ||
      !identical(fit$estimation_scope, "pathway_only") ||
      !identical(fit$scores_source, "unsmoothed_expression"))
    stop("fit must be an unsmoothed passage_fit_pathway() result; GP fits are not test inputs")
  if (!is.null(coords) && !is.null(rownames(coords)) &&
      !identical(rownames(coords), rownames(fit$scores)))
    stop("coordinate spot IDs must exactly match the fitted scores")
  if (!is.null(spatial_kernel) &&
      ((!is.null(rownames(spatial_kernel)) && !identical(rownames(spatial_kernel), rownames(fit$scores))) ||
       (!is.null(colnames(spatial_kernel)) && !identical(colnames(spatial_kernel), rownames(fit$scores)))))
    stop("spatial-kernel spot IDs must exactly match the fitted scores")
  ans <- passage_subspace_test(fit$scores, coords = coords, X = fit$X,
    spatial_kernel = spatial_kernel, kernel_range = kernel_range,
    K = fit$K, B = B, seed = seed, max_n = max_n)
  ans$estimation <- list(method = fit$method, scope = fit$estimation_scope,
    gene_names = fit$gene_names, K = fit$K, rank_rule = fit$rank_rule,
    explained_variance = fit$explained_variance, scores_source = fit$scores_source)
  # The inner calculation keeps every supplied score; report original gene-level PCA retention separately.
  ans$input_score_variance_retained <- ans$explained_variance
  ans$explained_variance <- fit$explained_variance
  ans$kernel_requirement <- "fixed or independently specified kernel; no same-data GP range selection"
  ans$calibration_status <- "Gaussian residual-space reference; real-data and FA performance require validation"
  ans
}
