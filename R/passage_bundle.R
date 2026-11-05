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
passage_sample_unique_sets <- function(pathway, strata, B = 9999L, seed = 1L,
                                       reference = c("complement", "universe")) {
  reference <- match.arg(reference)
  P <- unique(as.integer(pathway))
  if (!length(P) || anyNA(P) || any(P < 1 | P > length(strata))) stop("invalid pathway indices")
  if (anyNA(strata) || any(!nzchar(as.character(strata)))) stop("invalid matching strata")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  groups <- split(P, as.character(strata[P]))
  pools <- lapply(names(groups), function(s) {
    pool <- which(strata == s)
    if (reference == "complement") setdiff(pool, P) else pool
  })
  sizes <- lengths(groups)
  if (any(lengths(pools) < sizes)) stop("insufficient unique background genes in matching strata; coarsen prespecified bins")
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  lapply(seq_len(B), function(b) unlist(lapply(seq_along(pools), function(j)
    pools[[j]][sample.int(length(pools[[j]]), sizes[j], replace = FALSE)]), use.names = FALSE))
}

passage_analyze_bundle <- function(bundle, pathways = NULL, B = 9999L, seed = 1L,
                                   K = 6L, m = 20L, range_grid = NULL,
                                   gene_bins = NULL, min_size = 15L, max_size = 500L,
                                   verbose = TRUE, statistics = "score_z", combine = "none",
                                   rank_method = "fixed", variance_threshold = 0.90,
                                   max_K = NULL, min_K = 2L, factor_method = "pca",
                                   scale = FALSE, rank_grid = NULL,
                                   estimation = c("global", "pathway"),
                                   local_K = NULL, local_variance_threshold = 0.95,
                                   local_method = "pca", gp_aggregations = c("sum", "max"),
                                   gp_top_k = 2L, gp_kernel_method = "lowrank",
                                   gp_n_features = 64L, gp_ranges = c(0.05, 0.1, 0.2, 0.4),
                                   keep_matched_sets = FALSE,
                                   reference = c("complement", "universe")) {
  reference <- match.arg(reference)
  estimation <- match.arg(estimation)
  bundle <- passage_validate_bundle(bundle)
  Y <- bundle$Y
  keep <- apply(Y, 2L, stats::var) > 1e-12
  Y <- Y[, keep, drop = FALSE]
  if (is.null(pathways)) pathways <- bundle$pathway_sets
  if (!is.list(pathways) || !length(pathways) || is.null(names(pathways)) ||
      anyDuplicated(names(pathways)) || any(!nzchar(names(pathways))))
    stop("supply a named list of gene sets as pathways or bundle$pathway_sets")
  # Numeric definitions refer to the original bundle columns, before constants
  # are removed. Convert them to names so filtering cannot silently shift genes.
  pathways <- lapply(pathways, function(P) {
    if (is.numeric(P)) {
      if (anyNA(P) || any(!is.finite(P)) || any(P != as.integer(P)) ||
          any(P < 1L | P > ncol(bundle$Y))) stop("numeric pathway indices must lie within the original bundle gene universe")
      P <- colnames(bundle$Y)[P]
    }
    P
  })
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
  if (is.null(gene_bins)) gene_bins <- passage_make_gene_bins(Y,
    n_mean_bins = 3L, n_detect_bins = 2L, n_var_bins = 2L)
  if (!identical(as.character(gene_bins$gene), colnames(Y)))
    stop("gene_bins$gene must exactly match the nonconstant Y columns")
  # Check support for the entire declared family BEFORE fitting any models.
  matching_preflight <- lapply(idx, passage_matching_diagnostics,
    gene_bins = gene_bins, gene_names = colnames(Y), reference = reference)
  unsupported <- names(idx)[!vapply(matching_preflight, function(x) x$supported, logical(1))]
  if (length(unsupported)) stop("insufficient matching support for: ", paste(unsupported, collapse = ", "),
    "; use passage_matching_diagnostics() before selecting a prespecified binning scheme")
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  engine <- pre <- gene_scores <- NULL
  estimation_qc <- list()
  if (is.function(statistics)) {
    score_fun <- statistics
    estimation <- "custom"
  } else if (estimation == "pathway") {
    if (!is.null(rank_grid)) stop("rank_grid is for global estimation; prespecify local_K or local_variance_threshold")
    if (!identical(statistics, "score_z")) stop("use gp_aggregations for pathway estimation, or a custom statistic callback")
    score_fun <- passage_make_local_gp_statistic(Y, coords, X = bundle$X,
      method = local_method, K = local_K, variance_threshold = local_variance_threshold,
      ranges = gp_ranges, kernel_method = gp_kernel_method, n_features = gp_n_features,
      aggregations = gp_aggregations, top_k = gp_top_k, seed = seed)
    estimation_qc <- attr(score_fun, "configuration")
  } else {
    ranks <- if (is.null(rank_grid)) K else rank_grid
    if (!is.null(rank_grid) && (rank_method != "fixed" || anyNA(ranks) || !length(ranks) ||
        any(!is.finite(ranks)) || any(ranks < 1 | ranks != as.integer(ranks)) || anyDuplicated(ranks)))
      stop("rank_grid must contain distinct positive integer ranks with rank_method='fixed'")
    engines <- lapply(ranks, function(k) passage_fit_engine_pca(Y, coords, X = bundle$X,
      K = k, m = m, range_grid = range_grid, verbose = verbose,
      rank_method = rank_method, variance_threshold = variance_threshold,
      max_K = max_K, min_K = min_K, factor_method = factor_method, scale = scale))
    scorers <- lapply(engines, function(e) passage_make_statistic(e, Y, statistics = statistics,
      X = bundle$X, gene_bins = gene_bins))
    engine <- engines[[1L]]
    pre <- passage_h_precompute(engine, X = bundle$X)
    gene_scores <- passage_gene_score_z(engine, Y, precomp = pre)
    estimation_qc <- lapply(engines, function(e) list(K = e$K, rank_info = e$rank_info,
      factor_method = e$factor_method, theta = e$theta, scaled_genes = e$scale,
      design_columns = colnames(e$X), design_rank = qr(e$X)$rank))
    names(estimation_qc) <- paste0("K", vapply(engines, function(e) e$K, integer(1)))
    if (length(scorers) == 1L) score_fun <- scorers[[1L]] else {
      if (anyDuplicated(names(estimation_qc))) stop("rank_grid collapsed to duplicate fitted ranks; check max_K and residual degrees of freedom")
      score_fun <- function(P) {
        unlist(lapply(seq_along(scorers), function(i) {
          v <- scorers[[i]](P)
          stats::setNames(v, paste0(names(v), "_K", engines[[i]]$K))
        }), use.names = TRUE)
      }
      attr(score_fun, "configuration") <- list(estimation = "global", rank_grid = ranks,
        statistics = statistics, selection = "all ranks are separate components; combination calibration includes all of them")
    }
  }
  details <- vector("list", length(idx)); names(details) <- names(idx)
  rows <- vector("list", length(idx))
  test_rows <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    if (verbose) message("Competitive pathway ", j, "/", length(idx), ": ", names(idx)[j])
    P <- idx[[j]]
    test <- passage_competitive_test(P, score_fun, gene_bins, colnames(Y), B,
      if (is.null(seed)) NULL else seed + j, combine = combine,
      keep_matched_sets = keep_matched_sets, reference = reference)
    primary <- if (nrow(test$combinations)) test$combinations[1L, , drop = FALSE] else test$components[1L, , drop = FALSE]
    primary_name <- if (nrow(test$combinations)) paste0("combined_", primary$method) else primary$statistic
    obs <- primary$observed
    nul <- if (nrow(test$combinations)) test$combined_scores[-1L, 1L] else test$scores[-1L, 1L]
    p <- primary$p_competitive
    metrics <- if (!is.null(engine)) passage_pathway_covariance_metrics(engine, Y, P, precomp = pre)$summary else
      c(ePSV = NA_real_, cEPSV = NA_real_, mean_propSV_conditional = NA_real_, pc1_spatial_fraction = NA_real_, mean_abs_spatial_cor = NA_real_)
    driver <- if (!is.null(engine)) passage_sparse_topk_pathway_score_stat(engine, Y, P, precomp = pre,
      gene_scores = gene_scores, output = "details") else list(selected_genes = character())
    z <- if ("score_z" %in% colnames(test$scores)) test$scores[1L, "score_z"] else NA_real_
    rows[[j]] <- data.frame(pathway = names(idx)[j], n_genes = length(P), score_z = z,
      p_competitive = p, null_mean = mean(nul), null_sd = stats::sd(nul),
      score_excess = obs - mean(nul), mc_floor = 1 / (B + 1),
      ePSV = metrics[["ePSV"]], cEPSV = metrics[["cEPSV"]],
      mean_propSV = metrics[["mean_propSV_conditional"]],
      pc1_fraction = metrics[["pc1_spatial_fraction"]],
      coherence = metrics[["mean_abs_spatial_cor"]],
      effective_size = passage_effective_pathway_size(Y, P,
        precomp = if (is.null(pre)) list(qr = qr(bundle$X)) else pre),
      candidate_drivers = paste(driver$selected_genes, collapse = ";"),
      primary_test = primary_name, primary_statistic = obs,
      mc_lower = primary$mc_lower, mc_upper = primary$mc_upper, at_floor = primary$at_floor)
    tt <- test$components
    tt$test <- tt$statistic; tt$type <- "component"
    ct <- test$combinations
    if (nrow(ct)) {
      ct$test <- paste0("combined_", ct$method); ct$type <- "combination"
      common <- c("test", "type", "observed", "p_competitive", "exceedances", "B", "mc_lower", "mc_upper", "mc_floor", "at_floor")
      tt <- rbind(tt[, common], ct[, common])
    } else tt <- tt[, c("test", "type", "observed", "p_competitive", "exceedances", "B", "mc_lower", "mc_upper", "mc_floor", "at_floor")]
    tt$pathway <- names(idx)[j]
    test_rows[[j]] <- tt
    details[[j]] <- list(null = nul, metrics = metrics, candidate_drivers = driver, competitive = test)
  }
  summary <- do.call(rbind, rows)
  summary$q_BH <- stats::p.adjust(summary$p_competitive, "BH")
  summary$q_BY <- stats::p.adjust(summary$p_competitive, "BY")
  tests <- do.call(rbind, test_rows)
  tests$q_BH <- tests$q_BY <- NA_real_
  for (label in unique(tests$test)) {
    take <- tests$test == label
    tests$q_BH[take] <- stats::p.adjust(tests$p_competitive[take], "BH")
    tests$q_BY[take] <- stats::p.adjust(tests$p_competitive[take], "BY")
  }
  list(summary = summary, details = details, coverage = coverage,
       tests = tests, estimation = estimation_qc,
       configuration = list(version = "0.6.1", B = B, seed = seed, K = K, m = m,
         reference = reference,
         estimation = estimation, statistics = attr(score_fun, "configuration"), combine = combine,
         primary = "first requested combination, otherwise first statistic; prespecify this ordering",
         multiplicity = "summary across pathways for primary test; tests within each named test across pathways, not across statistic/combiner selection",
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
