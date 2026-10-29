# Validated, portable entry point for production input bundles.
passage_validate_bundle <- function(bundle, min_genes = 5L) {
  if (is.character(bundle) && length(bundle) == 1L) bundle <- readRDS(bundle)
  if (!is.list(bundle)) stop("bundle must be a list or an RDS path")
  Y <- passage_check_y(bundle$Y)
  if (nrow(Y) < 10L || ncol(Y) < min_genes) stop("bundle is too small")
  ids <- rownames(Y)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Y needs unique, nonempty spot row names")
  genes <- colnames(Y)
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes))
    stop("Y needs unique, nonempty gene names")
  coords <- passage_check_coords(bundle$coords, nrow(Y))
  if (!identical(rownames(coords), ids)) stop("coordinate spot IDs are not aligned")
  X <- passage_prepare_design(bundle$X, nrow(Y))
  if (!is.null(bundle$X) && !identical(rownames(X), ids)) stop("design spot IDs are not aligned")
  rownames(X) <- ids
  if (qr(X)$rank != ncol(X)) stop("design must have full column rank")
  if (nrow(Y) - ncol(X) < 5L) stop("fewer than five residual degrees of freedom")
  if (anyDuplicated(as.data.frame(coords))) stop("duplicate spatial coordinates")
  if (!is.null(bundle$cell_type_proportions)) {
    P <- as.matrix(bundle$cell_type_proportions)
    if (!identical(rownames(P), ids) || !all(is.finite(P)) || any(P < -1e-8) ||
        any(P > 1 + 1e-8) || any(abs(rowSums(P) - 1) > 1e-6))
      stop("cell_type_proportions must be aligned finite simplex fractions")
  }
  bundle$Y <- Y
  bundle$coords <- coords
  bundle$X <- X
  bundle
}

# A uniform draw within each prespecified stratum, without duplicate genes.
# Exclusion is strict: inadequate support is an error, never a fallback.
passage_sample_unique_sets <- function(pathway, strata, B = 9999L, seed = 1L) {
  P <- unique(as.integer(pathway))
  if (!length(P) || anyNA(P) || any(P < 1 | P > length(strata))) stop("invalid pathway indices")
  if (anyNA(strata) || any(!nzchar(as.character(strata)))) stop("invalid matching strata")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  groups <- split(P, as.character(strata[P]))
  pools <- lapply(names(groups), function(s) setdiff(which(strata == s), P))
  sizes <- lengths(groups)
  if (any(lengths(pools) < sizes)) stop("insufficient unique background genes in matching strata; coarsen prespecified bins")
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  lapply(seq_len(B), function(b) unlist(lapply(seq_along(pools), function(j)
    pools[[j]][sample.int(length(pools[[j]]), sizes[j], replace = FALSE)]), use.names = FALSE))
}

passage_analyze_bundle <- function(bundle, pathways = NULL, B = 9999L, seed = 1L,
                                   K = 6L, m = 20L, range_grid = NULL,
                                   gene_bins = NULL, min_size = 15L, max_size = 500L,
                                   verbose = TRUE) {
  bundle <- passage_validate_bundle(bundle)
  Y <- bundle$Y
  keep <- apply(Y, 2L, stats::var) > 1e-12
  Y <- Y[, keep, drop = FALSE]
  if (is.null(pathways)) pathways <- bundle$pathway_sets
  if (!is.list(pathways) || !length(pathways) || is.null(names(pathways)) ||
      anyDuplicated(names(pathways)) || any(!nzchar(names(pathways))))
    stop("supply a named list of gene sets as pathways or bundle$pathway_sets")
  idx <- lapply(pathways, passage_resolve_pathway, gene_names = colnames(Y))
  eligible <- lengths(idx) >= min_size & lengths(idx) <= max_size
  coverage <- data.frame(pathway = names(idx), n_input = lengths(pathways),
                         n_observed = lengths(idx), eligible = eligible)
  idx <- idx[eligible]
  if (!length(idx)) stop("no pathways pass the observed-gene size thresholds")
  coords <- bundle$coords
  # One scale for both axes preserves physical aspect ratio.
  coords <- sweep(coords, 2L, apply(coords, 2L, min), "-")
  coords <- coords / max(coords)
  set.seed(seed)
  engine <- passage_fit_engine_pca(Y, coords, X = bundle$X, K = K, m = m,
                                   range_grid = range_grid, verbose = verbose)
  pre <- passage_h_precompute(engine, X = bundle$X)
  if (is.null(gene_bins)) gene_bins <- passage_make_gene_bins(Y,
    n_mean_bins = 3L, n_detect_bins = 2L, n_var_bins = 2L)
  if (!identical(as.character(gene_bins$gene), colnames(Y)))
    stop("gene_bins$gene must exactly match the nonconstant Y columns")
  strata <- as.character(gene_bins$bin)
  if ("module" %in% names(gene_bins)) strata <- paste(strata, gene_bins$module, sep = "|")
  gene_scores <- passage_gene_score_z(engine, Y, precomp = pre)
  residuals <- passage_residualize_with_qr(Y, pre$qr)
  details <- vector("list", length(idx)); names(details) <- names(idx)
  rows <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    if (verbose) message("Competitive pathway ", j, "/", length(idx), ": ", names(idx)[j])
    P <- idx[[j]]
    matched <- passage_sample_unique_sets(P, strata, B, seed + j)
    all_scores <- passage_score_set_batch(engine, residuals, c(list(P), matched), pre)
    obs <- all_scores[1L]
    nul <- all_scores[-1L]
    if (!is.finite(obs) || any(!is.finite(nul))) stop("nonfinite observed/null statistic")
    p <- (1 + sum(nul >= obs)) / (B + 1)
    metrics <- passage_pathway_covariance_metrics(engine, Y, P, precomp = pre)$summary
    driver <- passage_sparse_topk_pathway_score_stat(engine, Y, P, precomp = pre,
      gene_scores = gene_scores, output = "details")
    rows[[j]] <- data.frame(pathway = names(idx)[j], n_genes = length(P), score_z = obs,
      p_competitive = p, null_mean = mean(nul), null_sd = stats::sd(nul),
      score_excess = obs - mean(nul), mc_floor = 1 / (B + 1),
      ePSV = metrics[["ePSV"]], cEPSV = metrics[["cEPSV"]],
      mean_propSV = metrics[["mean_propSV_conditional"]],
      pc1_fraction = metrics[["pc1_spatial_fraction"]],
      coherence = metrics[["mean_abs_spatial_cor"]],
      effective_size = passage_effective_pathway_size(Y, P, precomp = pre),
      candidate_drivers = paste(driver$selected_genes, collapse = ";"))
    details[[j]] <- list(null = nul, metrics = metrics, candidate_drivers = driver)
  }
  summary <- do.call(rbind, rows)
  summary$q_BH <- stats::p.adjust(summary$p_competitive, "BH")
  summary$q_BY <- stats::p.adjust(summary$p_competitive, "BY")
  list(summary = summary, details = details, coverage = coverage,
       configuration = list(version = "0.5.0", B = B, seed = seed, K = K, m = m,
         matching = if ("module" %in% names(gene_bins)) "unique_bin_and_module" else "unique_bin",
         calibration_status = "study-specific calibration required",
         driver_status = "candidate ranks; require spatial-block and held-out validation"),
       sample = bundle$sample, session_info = utils::sessionInfo())
}

# Algebraically identical to the equal-weight score_z, with shared residuals
# and sparse set-membership multiplication instead of repeatedly scanning Y.
passage_score_set_batch <- function(engine, residuals, sets, precomp, chunk_size = 64L) {
  out <- numeric(length(sets))
  weights <- passage_factor_weights(engine, "equal")
  for (first in seq.int(1L, length(sets), by = chunk_size)) {
    take <- first:min(length(sets), first + chunk_size - 1L)
    sub <- lapply(sets[take], unique)
    H <- Matrix::sparseMatrix(i = unlist(sub, use.names = FALSE),
      j = rep(seq_along(sub), lengths(sub)), x = 1,
      dims = c(ncol(residuals), length(sub)))
    Q <- mu <- vv <- numeric(length(sub))
    for (k in seq_len(engine$K)) {
      W <- Matrix::Diagonal(x = engine$A[, k]) %*% H
      Z <- residuals %*% W
      Z <- Z[engine$vecchia[[k]]$ord, , drop = FALSE]
      q <- as.numeric(Matrix::colSums(Z * (engine$K_score[[k]] %*% Z)))
      cvar <- as.numeric(crossprod(engine$A[, k]^2 * engine$D, H))
      Q <- Q + weights[k] * q
      mu <- mu + weights[k] * cvar * precomp$moments[[k]]$trace_mk
      vv <- vv + weights[k]^2 * 2 * cvar^2 * precomp$moments[[k]]$trace_mkmk
    }
    out[take] <- (Q - mu) / sqrt(pmax(vv, .Machine$double.eps))
  }
  out
}
