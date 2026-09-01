# Real-data-derived residual-null calibration for PASSAGE competitive tests.
#
# This script uses prepared breast Visium data, residualizes expression on the
# prepared covariate matrix, permutes residual rows with one common spot
# permutation across genes, and then computes matched competitive p-values.
#
# Usage:
#   Rscript scripts/calibrate_passage_competitive_real_residual_null.R \
#     --condition=H3 --n-reps=5 --B=99 --max-pathways=4 --pathway-source=hallmark \
#     --out-dir=results/passage_competitive_real_residual_null_quick

parse_args <- function(args) {
  cfg <- list(
    seed = 20260730L,
    n_reps = 5L,
    B = 99L,
    max_pathways = 4L,
    K = 3L,
    m = 8L,
    n_folds = 4L,
    modes = c("none", "spot_crossfit"),
    matching = c("expr_detect", "expr_detect_var", "expr_detect_var_spatial",
                 "expr_detect_var_spatial_corr",
                 "expr_detect_var_spatial_coherence",
                 "expr_detect_var_spatial_factor",
                 "expr_detect_var_spatial_corr_factor"),
    n_mean_bins = 10L,
    n_detect_bins = 10L,
    n_var_bins = 5L,
    n_spatial_bins = 5L,
    coherence = "signed",
    condition = "H1",
    score_statistic = "score_z",
    score_weight_scheme = "equal",
    mc_sampler = "independent",
    two_layer_gene_sampler = "module",
    n_modules = 30L,
    mcmc_burnin = 200L,
    mcmc_thin = 20L,
    importance_oversample = 20L,
    importance_temperature = 1,
    pathway_source = "hallmark",
    random_per_size = 1L,
    balance_rel_tol = 0.35,
    balance_abs_tol = 0.03,
    balance_max_tries = 50L,
    out_dir = file.path("results", "passage_competitive_real_residual_null_quick")
  )
  for (arg in args) {
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--n-reps=", arg)) cfg$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--B=", arg)) cfg$B <- as.integer(sub("^--B=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
    if (grepl("^--modes=", arg)) cfg$modes <- strsplit(sub("^--modes=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--matching=", arg)) cfg$matching <- strsplit(sub("^--matching=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--n-mean-bins=", arg)) cfg$n_mean_bins <- as.integer(sub("^--n-mean-bins=", "", arg))
    if (grepl("^--n-detect-bins=", arg)) cfg$n_detect_bins <- as.integer(sub("^--n-detect-bins=", "", arg))
    if (grepl("^--n-var-bins=", arg)) cfg$n_var_bins <- as.integer(sub("^--n-var-bins=", "", arg))
    if (grepl("^--n-spatial-bins=", arg)) cfg$n_spatial_bins <- as.integer(sub("^--n-spatial-bins=", "", arg))
    if (grepl("^--coherence=", arg)) cfg$coherence <- sub("^--coherence=", "", arg)
    if (grepl("^--condition=", arg)) cfg$condition <- toupper(sub("^--condition=", "", arg))
    if (grepl("^--score-statistic=", arg)) cfg$score_statistic <- sub("^--score-statistic=", "", arg)
    if (grepl("^--score-weight-scheme=", arg)) cfg$score_weight_scheme <- sub("^--score-weight-scheme=", "", arg)
    if (grepl("^--mc-sampler=", arg)) cfg$mc_sampler <- sub("^--mc-sampler=", "", arg)
    if (grepl("^--two-layer-gene-sampler=", arg)) cfg$two_layer_gene_sampler <- sub("^--two-layer-gene-sampler=", "", arg)
    if (grepl("^--n-modules=", arg)) cfg$n_modules <- as.integer(sub("^--n-modules=", "", arg))
    if (grepl("^--mcmc-burnin=", arg)) cfg$mcmc_burnin <- as.integer(sub("^--mcmc-burnin=", "", arg))
    if (grepl("^--mcmc-thin=", arg)) cfg$mcmc_thin <- as.integer(sub("^--mcmc-thin=", "", arg))
    if (grepl("^--importance-oversample=", arg)) cfg$importance_oversample <- as.integer(sub("^--importance-oversample=", "", arg))
    if (grepl("^--importance-temperature=", arg)) cfg$importance_temperature <- as.numeric(sub("^--importance-temperature=", "", arg))
    if (grepl("^--pathway-source=", arg)) cfg$pathway_source <- sub("^--pathway-source=", "", arg)
    if (grepl("^--random-per-size=", arg)) cfg$random_per_size <- as.integer(sub("^--random-per-size=", "", arg))
    if (grepl("^--balance-rel-tol=", arg)) cfg$balance_rel_tol <- as.numeric(sub("^--balance-rel-tol=", "", arg))
    if (grepl("^--balance-abs-tol=", arg)) cfg$balance_abs_tol <- as.numeric(sub("^--balance-abs-tol=", "", arg))
    if (grepl("^--balance-max-tries=", arg)) cfg$balance_max_tries <- as.integer(sub("^--balance-max-tries=", "", arg))
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
  }
  if (!cfg$pathway_source %in% c("hallmark", "random")) {
    stop("--pathway-source must be hallmark or random")
  }
  if (!is.finite(cfg$random_per_size) || cfg$random_per_size < 1L) {
    stop("--random-per-size must be a positive integer")
  }
  if (!cfg$condition %in% c("H1", "H3")) {
    stop("--condition must be H1 or H3")
  }
  if (!cfg$score_statistic %in% c("score_sparse_topk_z", "score_lse_z", "score_factor_hsic",
                                  "score_coherence", "score_activity_hotspot",
                                  "score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z",
                                  "score_pERSA_active_z", "score_robust_z", "score_z", "score_Q")) {
    stop("--score-statistic is not supported")
  }
  if (!cfg$mc_sampler %in% c("independent", "module", "mcmc", "importance", "knockoff", "two_layer")) {
    stop("--mc-sampler must be independent, module, mcmc, importance, knockoff, or two_layer")
  }
  if (!cfg$two_layer_gene_sampler %in% c("independent", "module", "mcmc", "importance", "knockoff")) {
    stop("--two-layer-gene-sampler must be independent, module, mcmc, importance, or knockoff")
  }
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

prepared_paths <- function(condition = "H1") {
  if (toupper(condition) == "H3") {
    return(c(
      Visium_FFPE_Human_Breast_Cancer =
        "results/passage_10x_hallmark_conditional_wu_perm999/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_conditional_prepared_data.rds",
      V1_Breast_Cancer_Block_A_Section_1 =
        "results/passage_10x_hallmark_conditional_wu_perm999/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_conditional_prepared_data.rds"
    ))
  }
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_perm999_fastcal/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_prepared_data.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_perm999_fastcal/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_prepared_data.rds"
  )
}

conditional_result_paths <- function() {
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_conditional_wu_perm999/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_conditional_result.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_conditional_wu_perm999/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_conditional_result.rds"
  )
}

condition_designs <- function(dat, dataset, cfg) {
  if (cfg$condition == "H1") {
    return(list(X_null = dat$X, X_score = dat$X))
  }
  fit <- readRDS(conditional_result_paths()[[dataset]])
  X3 <- cbind(as.matrix(dat$X), as.matrix(fit$Z_CT), as.matrix(fit$V_BG))
  list(X_null = X3, X_score = X3)
}

select_pathways <- function(pathways, max_pathways) {
  preferred <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA"
  )
  keep <- c(intersect(preferred, names(pathways)), setdiff(names(pathways), preferred))
  pathways[keep[seq_len(min(max_pathways, length(keep)))]]
}

make_random_pathways_like <- function(template_pathways, gene_names, seed, per_size = 1L) {
  set.seed(seed)
  sizes <- pmax(2L, pmin(lengths(template_pathways), length(gene_names)))
  out <- vector("list", length(sizes) * per_size)
  nn <- 0L
  for (ii in seq_along(sizes)) {
    for (rr in seq_len(per_size)) {
      nn <- nn + 1L
      out[[nn]] <- sample(gene_names, sizes[[ii]], replace = FALSE)
    }
  }
  names(out) <- unlist(lapply(seq_along(sizes), function(ii) {
    paste0("random_like_", ii, "_draw", seq_len(per_size), "_size", sizes[[ii]])
  }), use.names = FALSE)
  out
}

make_residual_null_y <- function(Y, X, seed) {
  set.seed(seed)
  Xd <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  qrx <- qr(Xd)
  R <- passage_residualize_with_qr(Y, qrx)
  Yhat <- Y - R
  Yhat + R[sample.int(nrow(R)), , drop = FALSE]
}

fit_engine_for_mode <- function(Y, coords, X, mode, rep_id, cfg) {
  range_grid <- passage_default_range_grid(coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
  fit_args <- list(
    K = cfg$K,
    m = cfg$m,
    range_grid = range_grid,
    kernel = "matern32",
    ordering = "none",
    verbose = FALSE
  )
  if (mode == "none") {
    return(passage_fit_engine_pca(
      Y = Y, coords = coords, X = X, K = cfg$K, m = cfg$m,
      range_grid = range_grid, kernel = "matern32", verbose = FALSE
    ))
  }
  if (mode == "spot_crossfit") {
    return(passage_fit_engine_crossfit_spots(
      Y = Y, coords = coords, X = X, n_folds = cfg$n_folds,
      seed = cfg$seed + rep_id, fit_args = fit_args, verbose = FALSE
    ))
  }
  if (mode == "pathway_holdout") {
    return(NULL)
  }
  stop("unsupported mode for competitive calibration: ", mode)
}

fit_pathway_holdout_engine <- function(Y, coords, X, genes, gene_names, cfg) {
  range_grid <- passage_default_range_grid(coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
  fit_args <- list(
    K = cfg$K,
    m = cfg$m,
    range_grid = range_grid,
    kernel = "matern32",
    ordering = "none",
    verbose = FALSE
  )
  passage_fit_engine_pathway_holdout(
    Y = Y,
    coords = coords,
    pathway = genes,
    X = X,
    gene_names = gene_names,
    fit_args = fit_args,
    verbose = FALSE
  )
}

make_bins_for_scheme <- function(Y, gene_names, gene_stats, scheme, cfg) {
  if (scheme == "expr_detect") {
    bins <- passage_make_gene_bins(
      Y, gene_names = gene_names,
      n_mean_bins = cfg$n_mean_bins,
      n_detect_bins = cfg$n_detect_bins,
      n_var_bins = 0L
    )
    bins$matching <- scheme
    return(bins)
  }
  if (scheme %in% c("expr_detect_var", "expr_detect_var_spatial",
                    "expr_detect_var_spatial_corr",
                    "expr_detect_var_spatial_coherence",
                    "expr_detect_var_spatial_factor",
                    "expr_detect_var_spatial_corr_factor")) {
    bins <- passage_make_gene_bins(
      Y, gene_names = gene_names,
      n_mean_bins = cfg$n_mean_bins,
      n_detect_bins = cfg$n_detect_bins,
      n_var_bins = cfg$n_var_bins
    )
    if (scheme %in% c("expr_detect_var_spatial", "expr_detect_var_spatial_corr",
                      "expr_detect_var_spatial_coherence",
                      "expr_detect_var_spatial_factor",
                      "expr_detect_var_spatial_corr_factor")) {
      bins$spatial_bin <- passage_rank_bins(gene_stats$propSV, cfg$n_spatial_bins)
      bins$bin <- paste(bins$bin, bins$spatial_bin, sep = ":")
    }
    bins$matching <- scheme
    return(bins)
  }
  stop("unknown matching scheme: ", scheme)
}

add_gene_modules <- function(bins, Y, pre, gene_stats, fast_ctx, cfg) {
  qrx <- if (!is.null(pre)) pre$qr else qr(passage_prepare_design(NULL, nrow(Y)))
  R <- passage_residualize_with_qr(Y, qrx)
  expr_pc <- tryCatch({
    s <- svd(scale(R, center = TRUE, scale = FALSE), nu = 0, nv = min(5L, ncol(R)))
    s$v
  }, error = function(e) matrix(0, nrow = ncol(Y), ncol = 1L))
  A <- as.matrix(fast_ctx$A)
  n_load <- min(5L, ncol(A))
  feat <- cbind(
    log_mean = log1p(bins$mean_expr),
    detection = bins$detection,
    log_var = log1p(pmax(bins$variance, 0)),
    propSV = gene_stats$propSV,
    expr_pc[, seq_len(min(5L, ncol(expr_pc))), drop = FALSE],
    A[, seq_len(n_load), drop = FALSE]
  )
  feat <- as.matrix(feat)
  feat[!is.finite(feat)] <- 0
  feat <- scale(feat)
  feat[!is.finite(feat)] <- 0
  k <- min(max(2L, as.integer(cfg$n_modules)), max(2L, nrow(feat) %/% 5L))
  set.seed(cfg$seed + 779L)
  km <- stats::kmeans(feat, centers = k, iter.max = 50L, nstart = 5L)
  bins$module <- paste0("module_", km$cluster)
  bins
}

mc_tail_p <- function(null, observed, weights = NULL, alternative = "greater") {
  ok <- is.finite(null)
  null <- null[ok]
  if (!is.null(weights)) {
    weights <- as.numeric(weights)[ok]
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(null))
  }
  if (length(null) == 0L) return(NA_real_)
  if (is.null(weights)) {
    return(switch(
      alternative,
      greater = (1 + sum(null >= observed)) / (1 + length(null)),
      less = (1 + sum(null <= observed)) / (1 + length(null)),
      stop("unsupported alternative")
    ))
  }
  tail <- switch(
    alternative,
    greater = sum(weights[null >= observed]) / sum(weights),
    less = sum(weights[null <= observed]) / sum(weights),
    stop("unsupported alternative")
  )
  (1 + length(null) * tail) / (1 + length(null))
}

score_pathway_competitive <- function(pathway_name, genes, gene_names, gene_stats,
                                      engine, Y, pre_score, fast_ctx, bins,
                                      expr_balance_fun, factor_balance_fun,
                                      dataset, mode, matching,
                                      rep_id, cfg) {
  P <- passage_resolve_pathway(genes, gene_names)
  t0 <- proc.time()[["elapsed"]]
  mean_score <- function(idx) mean(gene_stats$propSV[idx], na.rm = TRUE)
  score_value <- function(idx, output, Y_score = Y) {
    if (output == "robust_z") {
      gene_scores <- passage_gene_score_z(
        engine, Y_score, precomp = pre_score, gene_names = gene_names,
        weight_scheme = cfg$score_weight_scheme
      )
      return(passage_robust_pathway_score_stat(
        engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
        weight_scheme = cfg$score_weight_scheme, gene_scores = gene_scores
      ))
    }
    passage_pathway_score_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      weight_scheme = cfg$score_weight_scheme, output = output
    )
  }
  robust_gene_scores <- passage_gene_score_z(
    engine, Y, precomp = pre_score, gene_names = gene_names,
    weight_scheme = cfg$score_weight_scheme
  )
  robust_residuals <- passage_residualize_with_qr(Y, pre_score$qr)
  persa_options <- switch(
    cfg$score_statistic,
    score_pERSA_rank_z = list(reliability = "snr", standardization = "rank"),
    score_ePSA_rank_z = list(reliability = "none", standardization = "rank"),
    score_pERSA_active_z = list(reliability = "snr", standardization = "active"),
    list(reliability = "snr", standardization = "robust")
  )
  persa_gene_activity <- passage_gene_pERSA(
    engine, Y, precomp = pre_score, gene_names = gene_names,
    gene_bins = bins,
    reliability = persa_options$reliability,
    standardization = persa_options$standardization
  )
  lse_gene_activity <- passage_gene_lse(robust_gene_scores$score_z_gene, gene_bins = bins)
  spatial_features <- passage_residualize_with_qr(engine$V, pre_score$qr)
  robust_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    gene_scores <- if (use_cached) {
      robust_gene_scores
    } else {
      passage_gene_score_z(
        engine, Y_score, precomp = pre_score, gene_names = gene_names,
        weight_scheme = cfg$score_weight_scheme
      )
    }
    passage_robust_pathway_score_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      weight_scheme = cfg$score_weight_scheme, gene_scores = gene_scores,
      residual_matrix = if (use_cached) robust_residuals else NULL
    )
  }
  persa_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    gene_activity <- if (use_cached) {
      persa_gene_activity
    } else {
      passage_gene_pERSA(
        engine, Y_score, precomp = pre_score, gene_names = gene_names,
        gene_bins = bins,
        reliability = persa_options$reliability,
        standardization = persa_options$standardization
      )
    }
    passage_pERSA_pathway_score_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      gene_activity = gene_activity, gene_bins = bins,
      residual_matrix = if (use_cached) robust_residuals else NULL
    )
  }
  sparse_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    gene_scores <- if (use_cached) {
      robust_gene_scores
    } else {
      passage_gene_score_z(
        engine, Y_score, precomp = pre_score, gene_names = gene_names,
        weight_scheme = cfg$score_weight_scheme
      )
    }
    passage_sparse_topk_pathway_score_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      gene_scores = gene_scores,
      residual_matrix = if (use_cached) robust_residuals else NULL,
      weight_scheme = cfg$score_weight_scheme
    )
  }
  lse_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    lse <- if (use_cached) {
      lse_gene_activity
    } else {
      gs <- passage_gene_score_z(
        engine, Y_score, precomp = pre_score, gene_names = gene_names,
        weight_scheme = cfg$score_weight_scheme
      )
      passage_gene_lse(gs$score_z_gene, gene_bins = bins)
    }
    passage_gene_activity_pathway_score_stat(
      lse$lse_z, idx, Y_score, precomp = pre_score, gene_names = gene_names,
      residual_matrix = if (use_cached) robust_residuals else NULL
    )
  }
  factor_hsic_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    passage_pathway_factor_hsic_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      residual_matrix = if (use_cached) robust_residuals else NULL,
      spatial_features = spatial_features
    )
  }
  coherence_score <- function(idx) passage_fast_pc1_spatial_fraction(idx, fast_ctx)
  activity_hotspot_score <- function(idx, Y_score = Y) {
    use_cached <- identical(Y_score, Y)
    passage_pathway_activity_hotspot_stat(
      engine, Y_score, idx, precomp = pre_score, gene_names = gene_names,
      residual_matrix = if (use_cached) robust_residuals else NULL,
      weight_scheme = cfg$score_weight_scheme
    )
  }
  passage_score <- function(idx, Y_score = Y) {
    if (cfg$score_statistic == "score_Q") {
      score_value(idx, "Q", Y_score = Y_score)
    } else if (cfg$score_statistic == "score_sparse_topk_z") {
      sparse_score(idx, Y_score = Y_score)
    } else if (cfg$score_statistic == "score_lse_z") {
      lse_score(idx, Y_score = Y_score)
    } else if (cfg$score_statistic == "score_factor_hsic") {
      factor_hsic_score(idx, Y_score = Y_score)
    } else if (cfg$score_statistic == "score_coherence") {
      coherence_score(idx)
    } else if (cfg$score_statistic == "score_activity_hotspot") {
      activity_hotspot_score(idx, Y_score = Y_score)
    } else if (cfg$score_statistic %in% c("score_pERSA_z", "score_pERSA_rank_z",
                                          "score_ePSA_rank_z", "score_pERSA_active_z")) {
      persa_score(idx, Y_score = Y_score)
    } else if (cfg$score_statistic == "score_robust_z") {
      robust_score(idx, Y_score = Y_score)
    } else {
      score_value(idx, "z", Y_score = Y_score)
    }
  }
  cepsv_score <- function(idx) passage_fast_cEPSV(idx, fast_ctx)
  mean_obs <- mean_score(P)
  score_obs <- passage_score(P)
  score_q_obs <- NA_real_
  score_z_obs <- NA_real_
  score_robust_z_obs <- NA_real_
  score_pERSA_z_obs <- NA_real_
  score_sparse_topk_z_obs <- NA_real_
  score_lse_z_obs <- NA_real_
  score_factor_hsic_obs <- NA_real_
  score_coherence_obs <- NA_real_
  score_activity_hotspot_obs <- NA_real_
  if (cfg$score_statistic == "score_Q") {
    score_q_obs <- score_obs
  } else if (cfg$score_statistic == "score_z") {
    score_z_obs <- score_obs
  } else if (cfg$score_statistic == "score_robust_z") {
    score_robust_z_obs <- score_obs
  } else if (cfg$score_statistic %in% c("score_pERSA_z", "score_pERSA_rank_z",
                                        "score_ePSA_rank_z", "score_pERSA_active_z")) {
    score_pERSA_z_obs <- score_obs
  } else if (cfg$score_statistic == "score_sparse_topk_z") {
    score_sparse_topk_z_obs <- score_obs
  } else if (cfg$score_statistic == "score_lse_z") {
    score_lse_z_obs <- score_obs
  } else if (cfg$score_statistic == "score_factor_hsic") {
    score_factor_hsic_obs <- score_obs
  } else if (cfg$score_statistic == "score_coherence") {
    score_coherence_obs <- score_obs
  } else if (cfg$score_statistic == "score_activity_hotspot") {
    score_activity_hotspot_obs <- score_obs
  }
  cepsv_obs <- cepsv_score(P)
  balance_target <- NULL
  balance_tol <- NULL
  balance_max_tries <- 1L
  balance_fun <- NULL
  if (matching %in% c("expr_detect_var_spatial_corr", "expr_detect_var_spatial_coherence",
                      "expr_detect_var_spatial_factor",
                      "expr_detect_var_spatial_corr_factor")) {
    balance_fun <- switch(
      matching,
      expr_detect_var_spatial_corr = expr_balance_fun,
      expr_detect_var_spatial_coherence = expr_balance_fun,
      expr_detect_var_spatial_factor = factor_balance_fun,
      expr_detect_var_spatial_corr_factor = passage_combine_balance_functions(expr_balance_fun, factor_balance_fun)
    )
    balance_target <- balance_fun(P)
    if (matching == "expr_detect_var_spatial_corr") {
      balance_target <- balance_target["expr_mean_abs_cor"]
    } else if (matching == "expr_detect_var_spatial_factor") {
      balance_target <- balance_target[c("loading_pc1_fraction", "loading_mean_factor_coherence")]
    } else if (matching == "expr_detect_var_spatial_corr_factor") {
      balance_target <- balance_target[c("expr_mean_abs_cor",
                                         "loading_pc1_fraction",
                                         "loading_mean_factor_coherence")]
    }
    balance_tol <- pmax(abs(balance_target) * cfg$balance_rel_tol, cfg$balance_abs_tol)
    balance_max_tries <- cfg$balance_max_tries
  }
  seed_base <- cfg$seed + 100000L * rep_id + 1000L * match(mode, cfg$modes) +
    100L * match(matching, cfg$matching) + match(pathway_name, names(cfg$pathway_order))
  matched <- passage_sample_matched_indices(
    genes, bins, gene_names,
    B = cfg$B, seed = seed_base + 1L,
    sampler = if (cfg$mc_sampler == "two_layer") cfg$two_layer_gene_sampler else cfg$mc_sampler,
    mcmc_burnin = cfg$mcmc_burnin,
    mcmc_thin = cfg$mcmc_thin,
    importance_oversample = cfg$importance_oversample,
    importance_temperature = cfg$importance_temperature,
    balance_fun = balance_fun, balance_target = balance_target,
    balance_tol = balance_tol, balance_max_tries = balance_max_tries
  )
  mc_weights <- attr(matched, "weights")
  mean_null <- vapply(matched, mean_score, numeric(1))
  if (cfg$mc_sampler == "two_layer") {
    X_score <- pre_score$X
    score_null <- vapply(seq_along(matched), function(bb) {
      Yb <- make_residual_null_y(Y, X_score, seed_base + 1000000L + bb)
      passage_score(matched[[bb]], Y_score = Yb)
    }, numeric(1))
  } else {
    score_null <- vapply(matched, passage_score, numeric(1))
  }
  cepsv_null <- vapply(matched, cepsv_score, numeric(1))
  mean_null <- mean_null[is.finite(mean_null)]
  cepsv_null <- cepsv_null[is.finite(cepsv_null)]
  mean_p <- mc_tail_p(mean_null, mean_obs, weights = NULL)
  score_p <- mc_tail_p(score_null, score_obs, weights = mc_weights)
  cepsv_p <- mc_tail_p(cepsv_null, cepsv_obs, weights = NULL)
  mean_enrichment <- mean_obs / pmax(mean(mean_null), .Machine$double.eps)
  score_delta <- score_obs - mean(score_null[is.finite(score_null)])
  cepsv_enrichment <- cepsv_obs / pmax(mean(cepsv_null), .Machine$double.eps)
  balance <- attr(matched, "balance")
  balance_rate <- if (!is.null(balance)) balance$acceptance_rate else NA_real_
  balance_value <- function(nm) {
    if (!is.null(balance_target) && nm %in% names(balance_target)) balance_target[[nm]] else NA_real_
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  data.frame(
    dataset = dataset,
    replicate = rep_id,
    mode = mode,
    matching = matching,
    mc_sampler = cfg$mc_sampler,
    pathway = pathway_name,
    pathway_size = length(P),
    mean_propSV = mean_obs,
    competitive_mean_propSV_p = mean_p,
    competitive_mean_propSV_enrichment = mean_enrichment,
    score_Q = score_q_obs,
    score_z = score_z_obs,
    score_robust_z = score_robust_z_obs,
    score_pERSA_z = score_pERSA_z_obs,
    score_sparse_topk_z = score_sparse_topk_z_obs,
    score_lse_z = score_lse_z_obs,
    score_factor_hsic = score_factor_hsic_obs,
    score_coherence = score_coherence_obs,
    score_activity_hotspot = score_activity_hotspot_obs,
    competitive_score_p = score_p,
    competitive_score_delta = score_delta,
    cEPSV = cepsv_obs,
    competitive_cEPSV_p = cepsv_p,
    competitive_cEPSV_enrichment = cepsv_enrichment,
    balance_expr_mean_abs_cor = balance_value("expr_mean_abs_cor"),
    balance_expr_pc1_fraction = balance_value("expr_pc1_fraction"),
    balance_loading_pc1_fraction = balance_value("loading_pc1_fraction"),
    balance_loading_mean_factor_coherence = balance_value("loading_mean_factor_coherence"),
    mean_propSV_balance_acceptance = balance_rate,
    cEPSV_balance_acceptance = balance_rate,
    mc_weight_cv = if (!is.null(mc_weights)) stats::sd(mc_weights) / mean(mc_weights) else NA_real_,
    B = cfg$B,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_endpoint <- function(x, p_col) {
  groups <- unique(x[c("dataset", "mode", "matching", "mc_sampler")])
  rows <- list()
  ii <- 0L
  for (gg in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(x))
    for (cc in colnames(groups)) keep <- keep & x[[cc]] == groups[[cc]][gg]
    z <- x[keep, , drop = FALSE]
    for (alpha in c(0.10, 0.05, 0.01)) {
      p <- z[[p_col]][is.finite(z[[p_col]])]
      rate <- mean(p <= alpha)
      n <- length(p)
      se <- sqrt(rate * (1 - rate) / max(1L, n))
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        groups[gg, , drop = FALSE],
        endpoint = p_col,
        alpha = alpha,
        n = n,
        reject_rate = rate,
        mc_se = se,
        ci95_low = max(0, rate - 1.96 * se),
        ci95_high = min(1, rate + 1.96 * se),
        median_p = stats::median(p),
        min_p = min(p),
        median_elapsed_sec = stats::median(z$elapsed_sec),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

write_report <- function(cfg, summary) {
  md <- c(
    "# PASSAGE Competitive Real Residual Null Calibration",
    "",
    paste0("- Condition: ", cfg$condition),
    paste0("- Competitive score statistic: ", cfg$score_statistic),
    paste0("- Score weight scheme: ", cfg$score_weight_scheme),
    paste0("- Monte Carlo sampler: ", cfg$mc_sampler),
    paste0("- Gene modules: ", cfg$n_modules),
    paste0("- Reps per dataset: ", cfg$n_reps),
    paste0("- Pathways per replicate: ", cfg$max_pathways),
    paste0("- Matched gene-set permutations per pathway: ", cfg$B),
    paste0("- Modes: ", paste(cfg$modes, collapse = ", ")),
    paste0("- Matching schemes: ", paste(cfg$matching, collapse = ", ")),
    paste0("- Pathway source: ", cfg$pathway_source),
    paste0("- Random pathways per target size: ", cfg$random_per_size),
    paste0("- Set-balance relative tolerance: ", cfg$balance_rel_tol),
    paste0("- Set-balance absolute tolerance: ", cfg$balance_abs_tol),
    paste0("- Set-balance max tries: ", cfg$balance_max_tries),
    "",
    paste(colnames(summary), collapse = " | "),
    paste(rep("---", ncol(summary)), collapse = " | "),
    apply(summary, 1L, function(x) paste(x, collapse = " | "))
  )
  writeLines(md, file.path(cfg$out_dir, "summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
ii <- 0L
paths <- prepared_paths(cfg$condition)
for (dataset in names(paths)) {
  message("loading ", dataset)
  dat <- readRDS(paths[[dataset]])
  designs <- condition_designs(dat, dataset, cfg)
  gene_names <- colnames(dat$Y)
  template_pathways <- select_pathways(dat$pathways, cfg$max_pathways)
  if (cfg$pathway_source == "random") {
    pathways <- make_random_pathways_like(
      template_pathways, gene_names,
      cfg$seed + 10000L * match(dataset, names(paths)),
      per_size = cfg$random_per_size
    )
  } else {
    pathways <- template_pathways
  }
  cfg$pathway_order <- as.list(seq_along(pathways))
  names(cfg$pathway_order) <- names(pathways)
  for (rr in seq_len(cfg$n_reps)) {
    message("dataset=", dataset, " competitive residual-null rep=", rr)
    Y0 <- make_residual_null_y(dat$Y, designs$X_null, cfg$seed + 1000L * match(dataset, names(paths)) + rr)
    for (mode in cfg$modes) {
      message("  fitting mode=", mode)
      if (mode == "pathway_holdout") {
        for (pp in seq_along(pathways)) {
          message("    pathway_holdout pathway=", names(pathways)[pp])
          engine <- fit_pathway_holdout_engine(Y0, dat$coords, dat$X, pathways[[pp]], gene_names, cfg)
          pre <- passage_h_precompute(engine, X = passage_prepare_design(designs$X_score, nrow(Y0)))
          gene_stats <- passage_gene_spatial_stats(engine, Y0, precomp = pre, gene_names = gene_names)
          fast_ctx <- passage_competitive_fast_context(
            engine, Y0, precomp = pre, gene_names = gene_names, coherence = cfg$coherence
          )
          expr_balance_fun <- passage_make_expression_coherence_balance(
            Y0, precomp = pre, metrics = c("expr_mean_abs_cor", "expr_pc1_fraction")
          )
          factor_balance_fun <- passage_make_factor_coherence_balance(
            fast_ctx, metrics = c("loading_pc1_fraction", "loading_mean_factor_coherence")
          )
          for (matching in cfg$matching) {
            message("      matching=", matching)
            bins <- make_bins_for_scheme(Y0, gene_names, gene_stats, matching, cfg)
            bins <- add_gene_modules(bins, Y0, pre, gene_stats, fast_ctx, cfg)
            ii <- ii + 1L
            rows[[ii]] <- score_pathway_competitive(
              names(pathways)[pp], pathways[[pp]], gene_names, gene_stats,
              engine, Y0, pre, fast_ctx, bins, expr_balance_fun, factor_balance_fun,
              dataset, mode, matching, rr, cfg
            )
          }
        }
        next
      }
      engine <- fit_engine_for_mode(Y0, dat$coords, dat$X, mode, rr, cfg)
      pre <- passage_h_precompute(engine, X = passage_prepare_design(designs$X_score, nrow(Y0)))
      gene_stats <- passage_gene_spatial_stats(engine, Y0, precomp = pre, gene_names = gene_names)
      fast_ctx <- passage_competitive_fast_context(
        engine, Y0, precomp = pre, gene_names = gene_names, coherence = cfg$coherence
      )
      expr_balance_fun <- passage_make_expression_coherence_balance(
        Y0, precomp = pre, metrics = c("expr_mean_abs_cor", "expr_pc1_fraction")
      )
      factor_balance_fun <- passage_make_factor_coherence_balance(
        fast_ctx, metrics = c("loading_pc1_fraction", "loading_mean_factor_coherence")
      )
      for (matching in cfg$matching) {
        message("    matching=", matching)
        bins <- make_bins_for_scheme(Y0, gene_names, gene_stats, matching, cfg)
        bins <- add_gene_modules(bins, Y0, pre, gene_stats, fast_ctx, cfg)
        for (pp in seq_along(pathways)) {
          ii <- ii + 1L
          rows[[ii]] <- score_pathway_competitive(
            names(pathways)[pp], pathways[[pp]], gene_names, gene_stats,
            engine, Y0, pre, fast_ctx, bins, expr_balance_fun, factor_balance_fun,
            dataset, mode, matching, rr, cfg
          )
        }
      }
    }
  }
}

out <- do.call(rbind, rows)
out <- passage_empirical_competitive_calibration(out)
out$competitive_mean_propSV_fdr <- stats::p.adjust(out$competitive_mean_propSV_p, method = "BH")
out$competitive_score_fdr <- stats::p.adjust(out$competitive_score_p, method = "BH")
out$competitive_score_empirical_pooled_fdr <- stats::p.adjust(out$competitive_score_empirical_pooled_p, method = "BH")
out$competitive_score_empirical_leave_rep_fdr <- stats::p.adjust(out$competitive_score_empirical_leave_rep_p, method = "BH")
out$competitive_score_empirical_size_fdr <- stats::p.adjust(out$competitive_score_empirical_size_p, method = "BH")
out$competitive_score_empirical_leave_rep_size_fdr <- stats::p.adjust(out$competitive_score_empirical_leave_rep_size_p, method = "BH")
out$competitive_cEPSV_fdr <- stats::p.adjust(out$competitive_cEPSV_p, method = "BH")
write.csv(out, file.path(cfg$out_dir, "real_residual_null_competitive_pvalues.csv"), row.names = FALSE)

summary <- rbind(
  summarize_endpoint(out, "competitive_score_p"),
  summarize_endpoint(out, "competitive_score_empirical_pooled_p"),
  summarize_endpoint(out, "competitive_score_empirical_leave_rep_p"),
  summarize_endpoint(out, "competitive_score_empirical_size_p"),
  summarize_endpoint(out, "competitive_score_empirical_leave_rep_size_p"),
  summarize_endpoint(out, "competitive_mean_propSV_p"),
  summarize_endpoint(out, "competitive_cEPSV_p")
)
write.csv(summary, file.path(cfg$out_dir, "real_residual_null_competitive_summary.csv"), row.names = FALSE)
write_report(cfg, summary)
message("Wrote outputs to ", cfg$out_dir)
print(summary, row.names = FALSE)
