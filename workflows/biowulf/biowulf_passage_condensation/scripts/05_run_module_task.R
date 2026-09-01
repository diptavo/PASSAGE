#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
B <- if (length(args) >= 3L) as.integer(args[[3L]]) else 999L

source(file.path(root, "scripts", "passage_common.R"))
set.seed(20260817L + task_id * 811L)
tasks <- read.csv(file.path(root, "metadata", "module_task_manifest.csv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(tasks)) {
  message("Task outside manifest; exiting")
  quit(save = "no", status = 0)
}
task <- tasks[task_id, , drop = FALSE]
modules <- readRDS(file.path(root, "results", "pathway_modules", "pathway_module_gene_sets.rds"))
module_genes <- modules[[task$module_id[[1L]]]]
if (is.null(module_genes)) stop("Missing module gene set: ", task$module_id[[1L]])
pathways <- readRDS(file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds"))
background_genes <- sort(unique(unlist(pathways, use.names = FALSE)))

out_dir <- file.path(root, "results", "module_testing", task$cohort[[1L]], task$reference[[1L]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

obj <- readRDS(task$file[[1L]])
Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
Y <- Y[, intersect(colnames(Y), background_genes), drop = FALSE]
target_genes <- intersect(module_genes, colnames(Y))
if (length(target_genes) < 10L) stop("Too few module genes after overlap")
target_idx <- match(target_genes, colnames(Y))
ctx <- build_context(Y, as.matrix(obj$X), as.matrix(obj$coords))

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(ctx$moran8, 0)), sep = "_")

obs <- stat_all(target_idx, ctx)
null_mat <- matrix(NA_real_, nrow = B, ncol = length(obs), dimnames = list(NULL, names(obs)))
for (bb in seq_len(B)) {
  idx <- sample_matched(bins[target_idx], bins, length(target_idx), exclude = target_idx)
  null_mat[bb, ] <- stat_all(idx, ctx)
}

rows <- list()
rr <- 0L
for (st in names(obs)) {
  if (st != task$source_statistic[[1L]] && !(st %in% c("score_z", "score_z_robust", "CSPS", "GSPS", "HCPS"))) next
  tail <- gpd_tail_p(null_mat[, st], obs[[st]])
  rr <- rr + 1L
  rows[[rr]] <- data.frame(
    cohort = task$cohort[[1L]],
    sample = task$sample[[1L]],
    spatial_sample = task$spatial_sample[[1L]],
    reference = task$reference[[1L]],
    module_id = task$module_id[[1L]],
    source_statistic = task$source_statistic[[1L]],
    statistic = st,
    representative_pathway = task$representative_pathway[[1L]],
    n_member_pathways = task$n_member_pathways[[1L]],
      module_size = length(target_idx),
    observed = obs[[st]],
    null_mean = mean(null_mat[, st], na.rm = TRUE),
    null_sd = stats::sd(null_mat[, st], na.rm = TRUE),
    p_empirical = tail[["p_empirical"]],
    p_gpd = tail[["p_gpd"]],
    p_value = min(tail[["p_empirical"]], tail[["p_gpd"]], na.rm = TRUE),
    gpd_used = tail[["gpd_used"]],
    gpd_shape = tail[["gpd_shape"]],
    gpd_scale = tail[["gpd_scale"]],
    gpd_threshold = tail[["gpd_threshold"]],
    gpd_tail_n = tail[["gpd_tail_n"]],
    B = B,
      module_genes = paste(target_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
out <- do.call(rbind, rows)
prefix <- file.path(out_dir, sprintf("module_task_%06d", task_id))
write.csv(out, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(c(
  paste0("task_id=", task_id),
  paste0("module_id=", task$module_id[[1L]]),
  paste0("sample=", task$sample[[1L]]),
  paste0("B=", B),
  paste0("rows=", nrow(out)),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", paste0(prefix, ".csv"))
