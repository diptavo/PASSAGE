#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
max_pathways <- if (length(args) >= 2L) as.integer(args[[2L]]) else 80L
top_per_collection_stat <- if (length(args) >= 3L) as.integer(args[[3L]]) else 8L
q_cut <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 0.10

primary_stats <- c("score_z", "score_z_robust", "CSPS", "GSPS", "HCPS")
res_file <- file.path(root, "results", "pathway_testing_summary", "kidney_msigdb_pathway_all_results.csv")
if (!file.exists(res_file)) stop("Missing ", res_file)
res <- read.csv(res_file, stringsAsFactors = FALSE)
res <- res[res$statistic %in% primary_stats, , drop = FALSE]
task_tbl <- read.csv(file.path(root, "metadata", "pathway_task_manifest.csv"), stringsAsFactors = FALSE)

sample_rows <- list()
for (sample_id in sort(unique(res$sample))) {
  z0 <- res[res$sample == sample_id, , drop = FALSE]
  chosen <- character(0)
  for (coll in sort(unique(z0$collection))) {
    for (st in primary_stats) {
      z <- z0[z0$collection == coll & z0$statistic == st, , drop = FALSE]
      if (!nrow(z)) next
      z <- z[order(z$q_value_sample_collection_stat, z$p_value), , drop = FALSE]
      sig <- z$pathway[z$q_value_sample_collection_stat <= q_cut]
      top <- head(z$pathway, top_per_collection_stat)
      chosen <- unique(c(chosen, sig, top))
    }
  }
  p_by_path <- tapply(z0$p_value, z0$pathway, min, na.rm = TRUE)
  chosen <- chosen[chosen %in% names(p_by_path)]
  chosen <- chosen[order(p_by_path[chosen])]
  chosen <- head(unique(chosen), max_pathways)
  one <- z0[match(chosen, z0$pathway), c("sample", "sample_task_id", "reference"), drop = FALSE]
  sample_rows[[length(sample_rows) + 1L]] <- data.frame(
    driver_task_id = length(sample_rows) + 1L,
    sample = sample_id,
    sample_task_id = one$sample_task_id[[1L]],
    reference = one$reference[[1L]],
    file = unique(task_tbl$file[task_tbl$sample == sample_id])[1L],
    n_selected_pathways = length(chosen),
    selected_pathways = paste(chosen, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, sample_rows)
out_dir <- file.path(root, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(out_dir, "driver_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# Kidney PASSAGE Driver Manifest",
  "",
  paste0("- samples/references: ", nrow(out)),
  paste0("- max_pathways_per_sample: ", max_pathways),
  paste0("- top_per_collection_stat: ", top_per_collection_stat),
  paste0("- q_cut: ", q_cut),
  paste0("- primary_stats_for_pathway_selection: ", paste(primary_stats, collapse = ", "))
), file.path(out_dir, "driver_manifest.md"))
message("Wrote ", nrow(out), " driver tasks")
