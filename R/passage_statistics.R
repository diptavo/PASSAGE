# Statistic construction is independent of competitive null generation.
# All built-in statistics are oriented so that larger values are stronger.

.passage_statistic_names <- function() {
  c("score_z", "score_robust_z", "score_sparse_topk_z", "score_lse_z",
    "score_factor_hsic", "score_coherence", "score_activity_hotspot",
    "score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z",
    "score_pERSA_active_z", "score_Q", "mean_propSV_conditional",
    "cEPSV", "ePSV", "cwPVE_trace", "cwPVE_top", "mean_abs_spatial_cor")
}

# Recover expression in the same units as the fitted loadings and noise.
# Engine-specific scaling must not be applied to a cached residual matrix twice.
passage_engine_residuals <- function(engine, Y, precomp = NULL, X = NULL,
                                     indices = NULL) {
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (nrow(Y) != engine$N || ncol(Y) != engine$G)
    stop("Y dimensions must match the fitted engine")
  if (!is.null(engine$gene_names) && !identical(colnames(Y), engine$gene_names))
    stop("Y gene names and order must match the fitted engine")
  if (!is.null(rownames(engine$coords)) && !is.null(rownames(Y)) &&
      !identical(rownames(Y), rownames(engine$coords)))
    stop("Y spot names and order must match the fitted engine")
  if (is.null(indices)) indices <- seq_len(ncol(Y))
  if (!length(indices) || anyNA(indices) || any(indices != as.integer(indices)) ||
      any(indices < 1L | indices > ncol(Y))) stop("invalid gene indices")
  if (is.null(precomp)) {
    if (is.null(X)) X <- engine$X
    qrx <- qr(passage_prepare_design(X, nrow(Y), intercept = TRUE))
  } else {
    qrx <- precomp$qr
    if (is.null(qrx) || nrow(qrx$qr) != nrow(Y)) stop("precomp does not match Y")
  }
  R <- passage_residualize_with_qr(Y[, indices, drop = FALSE], qrx)
  if (!is.null(engine$gene_center))
    R <- sweep(R, 2L, engine$gene_center[indices], "-")
  if (!is.null(engine$gene_scale)) {
    sc <- engine$gene_scale[indices]
    if (any(!is.finite(sc) | sc <= 0)) stop("engine has invalid gene scales")
    R <- sweep(R, 2L, sc, "/")
  }
  R
}

passage_make_statistic <- function(engine, Y, statistics = "score_z",
                                    precomp = NULL, X = NULL, gene_bins = NULL,
                                    score_weight_scheme = "equal",
                                    covariance = "shrink", coherence = "signed",
                                    gene_names = NULL, ...) {
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  if (length(gene_names) != ncol(Y) || anyNA(gene_names) || anyDuplicated(gene_names))
    stop("gene_names must uniquely identify the columns of Y")
  statistics <- as.character(statistics)
  if (!length(statistics) || anyNA(statistics) || anyDuplicated(statistics) ||
      any(!statistics %in% .passage_statistic_names()))
    stop("statistics must be unique names from: ", paste(.passage_statistic_names(), collapse = ", "))
  score_weight_scheme <- match.arg(score_weight_scheme, c("equal", "var", "range", "invvar"))
  covariance <- match.arg(covariance, c("shrink", "sample", "diagonal"))
  coherence <- match.arg(coherence, c("signed", "absolute", "none"))
  dots <- list(...)
  allowed <- c("topk_grid", "winsor_q", "corr", "n_pc", "scale_floor_frac",
               "shrinkage", "ridge", "active_tail", "ridge_frac")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                      anyDuplicated(names(dots)) || any(!names(dots) %in% allowed)))
    stop("unsupported statistic option; allowed options: ", paste(allowed, collapse = ", "))
  opts <- function(keys) dots[intersect(names(dots), keys)]
  if (is.null(precomp)) precomp <- passage_h_precompute(engine, X = X)
  residuals <- passage_engine_residuals(engine, Y, precomp = precomp)
  uses_bins <- any(statistics %in% c("score_lse_z", "score_pERSA_z",
    "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z"))
  if (uses_bins && is.null(gene_bins)) gene_bins <- passage_make_gene_bins(Y, gene_names = gene_names)
  if (!is.null(gene_bins) && (!is.data.frame(gene_bins) || nrow(gene_bins) != ncol(Y) ||
       !"bin" %in% names(gene_bins) || anyNA(gene_bins$bin)))
    stop("gene_bins must be a data frame aligned with Y and contain bin")
  if (!is.null(gene_bins$gene) && !identical(as.character(gene_bins$gene), as.character(gene_names)))
    stop("gene_bins$gene must match gene_names in order")
  gene_scores <- if (any(statistics %in% c("score_robust_z", "score_sparse_topk_z", "score_lse_z")))
    passage_gene_score_z(engine, Y, precomp, gene_names, score_weight_scheme,
                         residual_matrix = residuals) else NULL
  lse <- if ("score_lse_z" %in% statistics)
    passage_gene_lse(gene_scores$score_z_gene, gene_bins = gene_bins) else NULL
  spatial_features <- if ("score_factor_hsic" %in% statistics)
    passage_residualize_with_qr(engine$V, precomp$qr) else NULL
  fast_ctx <- if ("score_coherence" %in% statistics)
    passage_competitive_fast_context(engine, Y, precomp = precomp, gene_names = gene_names) else NULL
  gene_activity <- list()
  for (st in intersect(statistics, c("score_pERSA_z", "score_pERSA_rank_z",
                                   "score_ePSA_rank_z", "score_pERSA_active_z"))) {
    standardization <- switch(st, score_pERSA_z = "robust", score_pERSA_active_z = "active", "rank")
    gene_activity[[st]] <- do.call(passage_gene_pERSA, c(list(
      engine = engine, Y = Y, precomp = precomp, gene_names = gene_names,
      gene_bins = gene_bins, reliability = if (st == "score_ePSA_rank_z") "none" else "snr",
      standardization = standardization, residual_matrix = residuals),
      opts(c("scale_floor_frac", "active_tail", "ridge_frac"))))
  }
  gene_spatial <- if ("mean_propSV_conditional" %in% statistics)
    passage_gene_spatial_stats(engine, Y, precomp = precomp, gene_names = gene_names,
                               residual_matrix = residuals) else NULL
  metric_names <- intersect(statistics, c("cEPSV", "ePSV", "cwPVE_trace", "cwPVE_top", "mean_abs_spatial_cor"))
  resolve <- function(idx) {
    if (is.numeric(idx) && (any(!is.finite(idx)) || any(idx != as.integer(idx)) ||
                            any(idx < 1L | idx > ncol(Y)))) stop("invalid gene indices")
    idx <- passage_resolve_pathway(idx, gene_names)
    if (!length(idx)) stop("pathway has no genes present in Y")
    idx
  }
  score_fun <- function(idx) {
    idx <- resolve(idx)
    # These helpers consume cached gene statistics rather than indexed engine
    # loadings. Local indexing avoids rescanning every gene for each null set.
    common <- list(engine = engine, Y = Y[, idx, drop = FALSE],
                   pathway = seq_along(idx), precomp = precomp,
                   gene_names = gene_names[idx],
                   residual_matrix = residuals[, idx, drop = FALSE])
    local_gene_scores <- if (!is.null(gene_scores)) gene_scores[idx, , drop = FALSE] else NULL
    joint <- if (any(c("score_z", "score_Q") %in% statistics))
      passage_score_from_residuals(engine, residuals[, idx, drop = FALSE], idx,
        engine$A[idx, , drop = FALSE], engine$D[idx], precomp,
        score_weight_scheme, FALSE, FALSE, integer(), gene_names)$joint[[score_weight_scheme]] else NULL
    metrics <- if (length(metric_names)) do.call(passage_pathway_covariance_metrics,
      c(list(engine = engine, Y = Y, pathway = idx, precomp = precomp,
             gene_names = gene_names, residual_matrix = residuals,
             covariance = covariance, coherence = coherence),
        opts(c("shrinkage", "scale_floor_frac", "ridge"))))$summary else NULL
    values <- vapply(statistics, function(st) switch(st,
      score_z = (joint$Q - joint$mean) / sqrt(pmax(joint$var, .Machine$double.eps)),
      score_Q = joint$Q,
      score_robust_z = do.call(passage_robust_pathway_score_stat,
        c(common, list(gene_scores = local_gene_scores, weight_scheme = score_weight_scheme), opts(c("winsor_q", "corr")))),
      score_sparse_topk_z = do.call(passage_sparse_topk_pathway_score_stat,
        c(common, list(gene_scores = local_gene_scores, weight_scheme = score_weight_scheme), opts("topk_grid"))),
      score_lse_z = do.call(passage_gene_activity_pathway_score_stat,
        c(common[setdiff(names(common), "engine")], list(activity = lse$lse_z[idx]), opts(c("winsor_q", "corr")))),
      score_factor_hsic = do.call(passage_pathway_factor_hsic_stat,
        c(common, list(spatial_features = spatial_features), opts("n_pc"))),
      score_coherence = passage_fast_pc1_spatial_fraction(idx, fast_ctx),
      score_activity_hotspot = do.call(passage_pathway_activity_hotspot_stat,
        c(common, list(weight_scheme = score_weight_scheme))),
      score_pERSA_z = , score_pERSA_rank_z = , score_ePSA_rank_z = , score_pERSA_active_z =
        do.call(passage_pERSA_pathway_score_stat,
          c(common, list(gene_activity = gene_activity[[st]][idx, , drop = FALSE]), opts(c("winsor_q", "corr")))),
      mean_propSV_conditional = mean(gene_spatial$propSV[idx]),
      metrics[[st]]), numeric(1))
    if (any(!is.finite(values))) stop("nonfinite pathway statistic: ",
      paste(names(values)[!is.finite(values)], collapse = ", "),
      "; inspect constant genes, pathway size, and fitted rank")
    values
  }
  if (identical(statistics, "score_z") && score_weight_scheme == "equal") {
    attr(score_fun, "batch") <- function(sets) {
      if (!is.list(sets) || !length(sets)) stop("sets must be a nonempty list")
      out <- matrix(passage_score_set_batch(engine, residuals, lapply(sets, resolve), precomp),
                    ncol = 1L, dimnames = list(names(sets), "score_z"))
      if (any(!is.finite(out))) stop("nonfinite batch statistic")
      out
    }
  }
  attr(score_fun, "configuration") <- list(statistics = statistics,
    score_weight_scheme = score_weight_scheme, covariance = covariance,
    coherence = coherence, statistic_options = dots, estimation = "global",
    K = engine$K, factor_method = engine$factor_method,
    scaled_genes = isTRUE(engine$scale), design_rank = precomp$qr$rank,
    alternative = "greater", calibration = "none; statistic only")
  score_fun
}
