#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
B <- if (length(args) >= 3L) as.integer(args[[3L]]) else 199L
min_genes <- if (length(args) >= 4L) as.integer(args[[4L]]) else 15L

source(file.path(root, "scripts", "passage_common.R"))
set.seed(20260817L + task_id * 1009L)
tasks <- read.csv(file.path(root, "metadata", "pathway_task_manifest.csv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(tasks)) {
  message("Task outside manifest; exiting")
  quit(save = "no", status = 0)
}
task <- tasks[task_id, , drop = FALSE]
pathways <- readRDS(file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds"))
pws <- strsplit(task$pathways[[1L]], ";", fixed = TRUE)[[1L]]
pws <- pws[pws %in% names(pathways)]

out_dir <- file.path(root, "results", "pathway_testing", task$cohort[[1L]], task$collection[[1L]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

obj <- readRDS(task$file[[1L]])
Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
pathway_genes <- sort(unique(unlist(pathways[pws], use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
if (ncol(Y) < 50L) stop("Too few genes after overlap")
ctx <- build_context(Y, as.matrix(obj$X), as.matrix(obj$coords))

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(ctx$moran8, 0)), sep = "_")

rows <- list()
rr <- 0L
for (pw in pws) {
  target_genes <- intersect(pathways[[pw]], colnames(Y))
  if (length(target_genes) < min_genes) next
  target_idx <- match(target_genes, colnames(Y))
  obs <- stat_all(target_idx, ctx)
  null_mat <- matrix(NA_real_, nrow = B, ncol = length(obs), dimnames = list(NULL, names(obs)))
  for (bb in seq_len(B)) {
    idx <- sample_matched(bins[target_idx], bins, length(target_idx), exclude = target_idx)
    null_mat[bb, ] <- stat_all(idx, ctx)
  }
  for (st in names(obs)) {
    tail <- gpd_tail_p(null_mat[, st], obs[[st]])
    rr <- rr + 1L
    rows[[rr]] <- data.frame(
      cohort = task$cohort[[1L]],
      sample = task$sample[[1L]],
      spatial_sample = task$spatial_sample[[1L]],
      reference = task$reference[[1L]],
      collection = task$collection[[1L]],
      pathway = pw,
      pathway_size = length(target_idx),
      statistic = st,
      observed = obs[[st]],
      null_mean = mean(null_mat[, st], na.rm = TRUE),
      null_sd = stats::sd(null_mat[, st], na.rm = TRUE),
      p_empirical = tail[["p_empirical"]],
      p_gpd = tail[["p_gpd"]],
      gpd_used = tail[["gpd_used"]],
      gpd_shape = tail[["gpd_shape"]],
      gpd_scale = tail[["gpd_scale"]],
      gpd_threshold = tail[["gpd_threshold"]],
      gpd_tail_n = tail[["gpd_tail_n"]],
      B = B,
      genes = paste(target_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (nrow(out)) {
  out$p_value <- pmin(out$p_empirical, out$p_gpd, na.rm = TRUE)
  out$q_value_sample_stat <- ave(out$p_value, out$sample, out$statistic, FUN = function(x) p.adjust(x, "BH"))
}
prefix <- file.path(out_dir, sprintf("pathway_task_%05d", task_id))
write.csv(out, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(c(
  paste0("task_id=", task_id),
  paste0("cohort=", task$cohort[[1L]]),
  paste0("sample=", task$sample[[1L]]),
  paste0("collection=", task$collection[[1L]]),
  paste0("rows=", nrow(out)),
  paste0("B=", B),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", paste0(prefix, ".csv"))
