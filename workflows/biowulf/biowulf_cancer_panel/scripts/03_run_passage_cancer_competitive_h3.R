#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
B <- if (length(args) >= 3L) as.integer(args[[3L]]) else 199L
cores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 2L

setwd(root)
suppressPackageStartupMessages({
  library(Matrix)
})
source(file.path(root, "scripts", "load_passage.R"))

manifest <- read.csv(file.path(root, "data", "passage_inputs", "passage_input_manifest.csv"),
                     stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  message("task_id ", task_id, " outside manifest rows=", nrow(manifest), "; exiting.")
  quit(save = "no", status = 0)
}
sample_id <- manifest$sample[[task_id]]
message("Competitive H3 cancer sample=", sample_id, " B=", B)
fit_path <- file.path(root, "results", "passage_cancer_h1_h3", sample_id, "passage_cancer_h1_h3_result.rds")
if (!file.exists(fit_path)) stop("Missing H1/H3 result: ", fit_path)
fit <- readRDS(fit_path)
Y <- fit$data$Y
engine <- fit$engine
pre_h3 <- fit$precomp$H3
pathways <- fit$data$pathways
gene_names <- colnames(Y)

gene_bins <- passage_make_gene_bins(Y, gene_names = gene_names, n_mean_bins = 10L, n_detect_bins = 10L, n_var_bins = 5L)
gene_stats <- passage_gene_spatial_stats(engine, Y, precomp = pre_h3, gene_names = gene_names)
gene_bins$spatial_bin <- passage_rank_bins(gene_stats$propSV, 5L)
gene_bins$bin <- paste(gene_bins$bin, gene_bins$spatial_bin, sep = ":")

fast_ctx <- passage_competitive_fast_context(engine, Y, precomp = pre_h3, gene_names = gene_names)
statistics <- c("score_z", "score_robust_z")

score_pathway <- function(ii, statistic) {
  pname <- names(pathways)[[ii]]
  genes <- pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    passage_conditional_competitive_test(
      engine = engine,
      Y = Y,
      pathway = genes,
      gene_bins = gene_bins,
      precomp = pre_h3,
      gene_names = gene_names,
      statistic = statistic,
      B = B,
      seed = 300000L + 10000L * task_id + 100L * ii + match(statistic, statistics),
      sampler = "module",
      score_weight_scheme = "equal"
    ),
    error = function(e) e
  )
  P <- passage_resolve_pathway(genes, gene_names)
  if (inherits(res, "error")) {
    return(data.frame(
      sample = sample_id, cancer = fit$cancer, spatial_sample = fit$spatial_sample,
      reference = fit$reference, pathway = pname, pathway_size = length(P),
      statistic = statistic, competitive_score = NA_real_, competitive_score_p = NA_real_,
      competitive_score_null_mean = NA_real_, competitive_score_null_sd = NA_real_,
      cEPSV = NA_real_, coherence_pc1 = NA_real_,
      error = conditionMessage(res), elapsed_sec = proc.time()[["elapsed"]] - t0,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    sample = sample_id,
    cancer = fit$cancer,
    spatial_sample = fit$spatial_sample,
    reference = fit$reference,
    pathway = pname,
    pathway_size = length(P),
    statistic = statistic,
    competitive_score = res$permutation$statistic,
    competitive_score_p = res$permutation$p,
    competitive_score_null_mean = res$permutation$null_mean,
    competitive_score_null_sd = res$permutation$null_sd,
    cEPSV = passage_fast_cEPSV(P, fast_ctx),
    coherence_pc1 = passage_fast_pc1_spatial_fraction(P, fast_ctx),
    error = "",
    elapsed_sec = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
}

jobs <- expand.grid(ii = seq_along(pathways), statistic = statistics, stringsAsFactors = FALSE)
rows <- parallel::mclapply(seq_len(nrow(jobs)), function(jj) {
  score_pathway(jobs$ii[[jj]], jobs$statistic[[jj]])
}, mc.cores = cores, mc.preschedule = FALSE)
tbl <- do.call(rbind, rows)
tbl$competitive_score_fdr <- ave(tbl$competitive_score_p, tbl$statistic,
                                 FUN = function(p) stats::p.adjust(p, method = "BH"))
tbl <- tbl[order(tbl$statistic, tbl$competitive_score_p), , drop = FALSE]

out_dir <- file.path(root, "results", "passage_cancer_competitive_h3", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tbl, file.path(out_dir, "passage_cancer_competitive_h3.csv"), row.names = FALSE)
saveRDS(list(sample = sample_id, cancer = fit$cancer, spatial_sample = fit$spatial_sample,
             reference = fit$reference, table = tbl, B = B, gene_bins = gene_bins),
        file.path(out_dir, "passage_cancer_competitive_h3.rds"))
message("Wrote ", out_dir)
