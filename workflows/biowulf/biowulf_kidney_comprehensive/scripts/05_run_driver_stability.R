#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_boot <- if (length(args) >= 3L) as.integer(args[[3L]]) else 30L
n_null <- if (length(args) >= 4L) as.integer(args[[4L]]) else 9L
top_k <- if (length(args) >= 5L) as.integer(args[[5L]]) else 10L

source(file.path(root, "scripts", "passage_kidney_common.R"))
stat_names <- c("CSPS", "GSPS", "MMP", "CSV", "HCPS", "OTSAS", "score_z", "score_z_robust")
set.seed(20260816L + task_id * 997L)

manifest <- read.csv(file.path(root, "metadata", "driver_manifest.csv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  message("Task outside manifest; exiting")
  quit(save = "no", status = 0)
}
task <- manifest[task_id, , drop = FALSE]
pathways <- readRDS(file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds"))
pathways_to_run <- strsplit(task$selected_pathways[[1L]], ";", fixed = TRUE)[[1L]]
pathways_to_run <- pathways_to_run[nzchar(pathways_to_run) & pathways_to_run %in% names(pathways)]

out_dir <- file.path(root, "results", "driver_stability", task$reference[[1L]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

obj <- readRDS(task$file[[1L]])
Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
pathway_genes <- sort(unique(unlist(pathways[pathways_to_run], use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
X <- as.matrix(obj$X)
coords <- as.matrix(obj$coords)
full_ctx <- build_context(Y, X, coords)

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(full_ctx$moran8, 0)), sep = "_")

boot_ctx <- vector("list", n_boot)
for (bb in seq_len(n_boot)) boot_ctx[[bb]] <- build_context(Y, X, coords, keep = spatial_block_keep(coords))

driver_rows <- list()
summary_rows <- list()
null_rows <- list()
dr <- sr <- nr <- 0L
for (pw in pathways_to_run) {
  target_genes <- intersect(pathways[[pw]], colnames(Y))
  if (length(target_genes) < 10L) next
  target_idx <- match(target_genes, colnames(Y))
  null_sets <- vector("list", n_null)
  for (nn in seq_len(n_null)) null_sets[[nn]] <- sample_matched(bins[target_idx], bins, length(target_idx), exclude = target_idx)

  for (stat in stat_names) {
    selected_sets <- vector("list", n_boot)
    for (bb in seq_len(n_boot)) {
      scores <- driver_scores_for_idx(target_idx, boot_ctx[[bb]], stat)
      sel <- asset_select(scores, names(scores), top_k = top_k)
      selected_sets[[bb]] <- sel$gene
      sel$bootstrap <- bb
      sel$sample <- task$sample[[1L]]
      sel$reference <- task$reference[[1L]]
      sel$pathway <- pw
      sel$pathway_size <- length(target_idx)
      sel$statistic <- stat
      dr <- dr + 1L
      driver_rows[[dr]] <- sel[, c("sample", "reference", "pathway", "pathway_size", "statistic", "bootstrap", "gene", "rank", "driver_weight", "raw_score")]
    }
    freq <- sort(table(unlist(selected_sets, use.names = FALSE)) / n_boot, decreasing = TRUE)
    sr <- sr + 1L
    summary_rows[[sr]] <- data.frame(
      sample = task$sample[[1L]],
      reference = task$reference[[1L]],
      pathway = pw,
      pathway_size = length(target_idx),
      statistic = stat,
      n_boot = n_boot,
      top_k = top_k,
      n_unique_selected = length(freq),
      max_selection_frequency = if (length(freq)) as.numeric(freq[[1L]]) else NA_real_,
      n_genes_freq_ge_050 = sum(freq >= 0.50),
      n_genes_freq_ge_025 = sum(freq >= 0.25),
      mean_topk_jaccard = jaccard_mean(selected_sets),
      top_genes = paste(names(head(freq, 20L)), collapse = ";"),
      top_gene_frequencies = paste(sprintf("%.3f", as.numeric(head(freq, 20L))), collapse = ";"),
      stringsAsFactors = FALSE
    )
    for (nn in seq_len(n_null)) {
      nsets <- vector("list", n_boot)
      for (bb in seq_len(n_boot)) {
        nscores <- driver_scores_for_idx(null_sets[[nn]], boot_ctx[[bb]], stat)
        nsets[[bb]] <- asset_select(nscores, names(nscores), top_k = top_k)$gene
      }
      nf <- sort(table(unlist(nsets, use.names = FALSE)) / n_boot, decreasing = TRUE)
      nr <- nr + 1L
      null_rows[[nr]] <- data.frame(
        sample = task$sample[[1L]],
        reference = task$reference[[1L]],
        pathway = pw,
        statistic = stat,
        null_id = nn,
        null_max_selection_frequency = if (length(nf)) as.numeric(nf[[1L]]) else NA_real_,
        null_mean_topk_jaccard = jaccard_mean(nsets),
        stringsAsFactors = FALSE
      )
    }
  }
}

drivers <- if (length(driver_rows)) do.call(rbind, driver_rows) else data.frame()
summ <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
nulls <- if (length(null_rows)) do.call(rbind, null_rows) else data.frame()
if (nrow(summ) && nrow(nulls)) {
  key <- paste(summ$sample, summ$pathway, summ$statistic, sep = "\r")
  nkey <- paste(nulls$sample, nulls$pathway, nulls$statistic, sep = "\r")
  summ$null_p_max_frequency <- vapply(seq_len(nrow(summ)), function(i) {
    z <- nulls$null_max_selection_frequency[nkey == key[[i]]]
    (1 + sum(z >= summ$max_selection_frequency[[i]], na.rm = TRUE)) / (1 + length(z))
  }, numeric(1))
  summ$null_p_jaccard <- vapply(seq_len(nrow(summ)), function(i) {
    z <- nulls$null_mean_topk_jaccard[nkey == key[[i]]]
    (1 + sum(z >= summ$mean_topk_jaccard[[i]], na.rm = TRUE)) / (1 + length(z))
  }, numeric(1))
}

prefix <- file.path(out_dir, sprintf("driver_stability_task_%02d", task_id))
write.csv(drivers, paste0(prefix, "_drivers.csv"), row.names = FALSE)
write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)
write.csv(nulls, paste0(prefix, "_nulls.csv"), row.names = FALSE)
writeLines(c(
  paste0("task_id=", task_id),
  paste0("sample=", task$sample[[1L]]),
  paste0("pathways=", length(pathways_to_run)),
  paste0("statistics=", paste(stat_names, collapse = ",")),
  paste0("n_boot=", n_boot),
  paste0("n_null=", n_null),
  paste0("top_k=", top_k),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", prefix)
