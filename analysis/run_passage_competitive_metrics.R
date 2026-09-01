# Compute covariance-aware PASSAGE metrics and matched competitive tests
# from an existing PASSAGE conditional result directory.

suppressPackageStartupMessages({
  library(parallel)
})

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

parse_args <- function(args) {
  cfg <- list(
    result_dir = NULL,
    condition = "H3",
    B = 999L,
    cores = max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L)),
    seed = 20260524L,
    out = NULL,
    n_mean_bins = 10L,
    n_detect_bins = 10L,
    n_var_bins = 0L,
    covariance = "shrink",
    coherence = "signed",
    max_pathways = NA_integer_
  )
  for (arg in args) {
    if (grepl("^--result-dir=", arg)) cfg$result_dir <- sub("^--result-dir=", "", arg)
    if (grepl("^--condition=", arg)) cfg$condition <- sub("^--condition=", "", arg)
    if (grepl("^--B=", arg)) cfg$B <- as.integer(sub("^--B=", "", arg))
    if (grepl("^--cores=", arg)) cfg$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out=", arg)) cfg$out <- sub("^--out=", "", arg)
    if (grepl("^--n-mean-bins=", arg)) cfg$n_mean_bins <- as.integer(sub("^--n-mean-bins=", "", arg))
    if (grepl("^--n-detect-bins=", arg)) cfg$n_detect_bins <- as.integer(sub("^--n-detect-bins=", "", arg))
    if (grepl("^--n-var-bins=", arg)) cfg$n_var_bins <- as.integer(sub("^--n-var-bins=", "", arg))
    if (grepl("^--covariance=", arg)) cfg$covariance <- sub("^--covariance=", "", arg)
    if (grepl("^--coherence=", arg)) cfg$coherence <- sub("^--coherence=", "", arg)
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
  }
  if (is.null(cfg$result_dir)) stop("--result-dir is required")
  if (is.null(cfg$out)) {
    cfg$out <- file.path(cfg$result_dir, paste0("passage_competitive_metrics_", cfg$condition, ".csv"))
  }
  cfg
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
fit_path <- file.path(cfg$result_dir, "passage_msigdb_conditional_result.rds")
dat_path <- file.path(cfg$result_dir, "passage_msigdb_conditional_prepared_data.rds")
if (!file.exists(fit_path)) stop("Missing result RDS: ", fit_path)
if (!file.exists(dat_path)) stop("Missing prepared data RDS: ", dat_path)

fit <- readRDS(fit_path)
dat <- readRDS(dat_path)
if (is.null(fit$precomp[[cfg$condition]])) stop("Condition not found in fit$precomp: ", cfg$condition)
precomp <- fit$precomp[[cfg$condition]]
gene_names <- colnames(dat$Y)
gene_bins <- passage_make_gene_bins(
  dat$Y,
  gene_names = gene_names,
  n_mean_bins = cfg$n_mean_bins,
  n_detect_bins = cfg$n_detect_bins,
  n_var_bins = cfg$n_var_bins
)
gene_stats <- passage_gene_spatial_stats(fit$engine, dat$Y, precomp = precomp, gene_names = gene_names)
fast_ctx <- passage_competitive_fast_context(
  fit$engine, dat$Y, precomp = precomp, gene_names = gene_names,
  coherence = cfg$coherence
)
pathways <- dat$pathways
if (is.finite(cfg$max_pathways)) {
  pathways <- pathways[seq_len(min(length(pathways), cfg$max_pathways))]
}

message("Computing competitive PASSAGE metrics for ", length(pathways), " pathways")
message("Condition: ", cfg$condition, "; B=", cfg$B, "; cores=", cfg$cores)

score_one <- function(ii) {
  pname <- names(pathways)[ii]
  genes <- pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch({
    metrics <- passage_pathway_covariance_metrics(
      fit$engine, dat$Y, genes, precomp = precomp, gene_names = gene_names,
      covariance = cfg$covariance, coherence = cfg$coherence
    )$summary
    z <- passage_competitive_gene_stat_z(genes, gene_stats$propSV, gene_bins, gene_names)
    score_mean <- function(idx) mean(gene_stats$propSV[idx], na.rm = TRUE)
    mean_perm <- passage_competitive_permutation_test(
      genes, score_mean, gene_bins, gene_names,
      B = cfg$B, seed = cfg$seed + 1000L * ii + 1L,
      observed = z$statistic
    )
    cepsv_score <- function(idx) passage_fast_cEPSV(idx, fast_ctx)
    cepsv <- passage_competitive_permutation_test(
      genes, cepsv_score, gene_bins, gene_names,
      B = cfg$B, seed = cfg$seed + 1000L * ii + 2L,
      observed = metrics[["cEPSV"]]
    )
    coh_score <- function(idx) passage_fast_pc1_spatial_fraction(idx, fast_ctx)
    coh <- passage_competitive_permutation_test(
      genes, coh_score, gene_bins, gene_names,
      B = cfg$B, seed = cfg$seed + 1000L * ii + 3L,
      observed = metrics[["pc1_spatial_fraction"]]
    )
    row <- data.frame(
      pathway = pname,
      pathway_size = length(passage_resolve_pathway(genes, gene_names)),
      condition = cfg$condition,
      mean_propSV_conditional = metrics[["mean_propSV_conditional"]],
      cwPVE_trace = metrics[["cwPVE_trace"]],
      cwPVE_top = metrics[["cwPVE_top"]],
      cwPVE_mean = metrics[["cwPVE_mean"]],
      ePSV = metrics[["ePSV"]],
      cEPSV = metrics[["cEPSV"]],
      spatial_eff_rank = metrics[["spatial_eff_rank"]],
      residual_eff_rank = metrics[["residual_eff_rank"]],
      residual_mean_cor = metrics[["residual_mean_cor"]],
      residual_eff_genes_vif = metrics[["residual_eff_genes_vif"]],
      pc1_spatial_fraction = metrics[["pc1_spatial_fraction"]],
      mean_abs_spatial_cor = metrics[["mean_abs_spatial_cor"]],
      competitive_mean_propSV_z = z$z,
      competitive_mean_propSV_p_analytic = z$p,
      competitive_mean_propSV_expected = z$expected,
      competitive_mean_propSV_p_perm = mean_perm$p,
      competitive_mean_propSV_enrichment = mean_perm$enrichment,
      competitive_cEPSV_p_perm = cepsv$p,
      competitive_cEPSV_enrichment = cepsv$enrichment,
      coherence_pc1_p_perm = coh$p,
      coherence_pc1_enrichment = coh$enrichment,
      elapsed_sec = proc.time()[["elapsed"]] - t0,
      error_message = "",
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "COMPETITIVE_DONE\t%s\tcEPSV=%.4g\tcompP=%.4g\tcohP=%.4g\telapsed=%.1fs",
      pname, row$cEPSV, row$competitive_cEPSV_p_perm, row$coherence_pc1_p_perm, row$elapsed_sec
    ))
    row
  }, error = function(e) {
    message("COMPETITIVE_ERROR\t", pname, "\t", conditionMessage(e))
    data.frame(
      pathway = pname,
      pathway_size = length(passage_resolve_pathway(genes, gene_names)),
      condition = cfg$condition,
      mean_propSV_conditional = NA_real_,
      cwPVE_trace = NA_real_,
      cwPVE_top = NA_real_,
      cwPVE_mean = NA_real_,
      ePSV = NA_real_,
      cEPSV = NA_real_,
      spatial_eff_rank = NA_real_,
      residual_eff_rank = NA_real_,
      residual_mean_cor = NA_real_,
      residual_eff_genes_vif = NA_real_,
      pc1_spatial_fraction = NA_real_,
      mean_abs_spatial_cor = NA_real_,
      competitive_mean_propSV_z = NA_real_,
      competitive_mean_propSV_p_analytic = NA_real_,
      competitive_mean_propSV_expected = NA_real_,
      competitive_mean_propSV_p_perm = NA_real_,
      competitive_mean_propSV_enrichment = NA_real_,
      competitive_cEPSV_p_perm = NA_real_,
      competitive_cEPSV_enrichment = NA_real_,
      coherence_pc1_p_perm = NA_real_,
      coherence_pc1_enrichment = NA_real_,
      elapsed_sec = proc.time()[["elapsed"]] - t0,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  out
}

rows <- parallel::mclapply(seq_along(pathways), score_one, mc.cores = cfg$cores, mc.preschedule = FALSE)
tbl <- do.call(rbind, rows)
tbl$competitive_mean_propSV_fdr_analytic <- stats::p.adjust(tbl$competitive_mean_propSV_p_analytic, method = "BH")
tbl$competitive_mean_propSV_fdr_perm <- stats::p.adjust(tbl$competitive_mean_propSV_p_perm, method = "BH")
tbl$competitive_cEPSV_fdr_perm <- stats::p.adjust(tbl$competitive_cEPSV_p_perm, method = "BH")
tbl$coherence_pc1_fdr_perm <- stats::p.adjust(tbl$coherence_pc1_p_perm, method = "BH")
tbl <- tbl[order(tbl$competitive_cEPSV_p_perm, tbl$coherence_pc1_p_perm, na.last = TRUE), , drop = FALSE]
rownames(tbl) <- NULL
write.csv(tbl, cfg$out, row.names = FALSE)
message("Wrote ", cfg$out)
