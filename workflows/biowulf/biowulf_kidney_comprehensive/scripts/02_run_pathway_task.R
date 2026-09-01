#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
B <- if (length(args) >= 3L) as.integer(args[[3L]]) else 199L

source(file.path(root, "scripts", "passage_kidney_common.R"))
set.seed(20260816L + task_id * 1009L)

task_file <- file.path(root, "metadata", "pathway_task_manifest.csv")
pathway_file <- file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds")
task_tbl <- read.csv(task_file, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(task_tbl)) {
  message("Task outside manifest; exiting")
  quit(save = "no", status = 0)
}
task <- task_tbl[task_id, , drop = FALSE]
pathways <- readRDS(pathway_file)
pws <- strsplit(task$pathways[[1L]], ";", fixed = TRUE)[[1L]]
pws <- pws[pws %in% names(pathways)]

out_dir <- file.path(root, "results", "pathway_testing", task$collection[[1L]])
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
  if (length(target_genes) < 10L) next
  target_idx <- match(target_genes, colnames(Y))
  size <- length(target_idx)
  obs <- stat_all(target_idx, ctx)
  null_mat <- matrix(NA_real_, nrow = B, ncol = length(obs), dimnames = list(NULL, names(obs)))
  target_bins <- bins[target_idx]
  for (bb in seq_len(B)) {
    idx <- sample_matched(target_bins, bins, size, exclude = target_idx)
    null_mat[bb, ] <- stat_all(idx, ctx)
  }
  for (st in names(obs)) {
    rr <- rr + 1L
    rows[[rr]] <- data.frame(
      sample = task$sample[[1L]],
      sample_task_id = task$sample_task_id[[1L]],
      reference = sub("^kidney__KC[0-9]+__", "", task$sample[[1L]]),
      collection = task$collection[[1L]],
      pathway = pw,
      pathway_size = size,
      statistic = st,
      observed = obs[[st]],
      null_mean = mean(null_mat[, st], na.rm = TRUE),
      null_sd = stats::sd(null_mat[, st], na.rm = TRUE),
      p_value = (1 + sum(null_mat[, st] >= obs[[st]], na.rm = TRUE)) / (B + 1),
      B = B,
      genes = paste(target_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

out <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (nrow(out)) {
  out$q_value_sample_stat <- ave(out$p_value, out$sample, out$statistic, FUN = function(x) p.adjust(x, "BH"))
  out <- out[order(out$statistic, out$q_value_sample_stat, out$p_value), , drop = FALSE]
}
prefix <- file.path(out_dir, sprintf("pathway_task_%04d", task_id))
write.csv(out, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(c(
  paste0("task_id=", task_id),
  paste0("sample=", task$sample[[1L]]),
  paste0("collection=", task$collection[[1L]]),
  paste0("pathways_requested=", length(pws)),
  paste0("rows=", nrow(out)),
  paste0("B=", B),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", paste0(prefix, ".csv"))
