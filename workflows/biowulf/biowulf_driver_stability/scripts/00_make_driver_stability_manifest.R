#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
pathway_root <- if (length(args) >= 2L) args[[2L]] else "/data/DCEG_Dutta/PASSAGE_pathway_new_stats_20260814"
max_pathways <- if (length(args) >= 3L) as.integer(args[[3L]]) else 12L

cohort_roots <- c(
  breast = "/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803",
  kidney = "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803",
  dlpfc = "/data/Dutta_lab/SPATH/PASSAGE_spatialDLPFC_20260802"
)

path_file <- file.path(pathway_root, "results", "pathway_new_stats_summary", "pathway_new_stats_all_results.csv")
if (!file.exists(path_file)) stop("Missing pathway result table: ", path_file)
res <- read.csv(path_file, stringsAsFactors = FALSE)

sample_rows <- list()
for (cohort in names(cohort_roots)) {
  manifest <- read.csv(file.path(cohort_roots[[cohort]], "data", "passage_inputs", "passage_input_manifest.csv"),
                       stringsAsFactors = FALSE)
  if (cohort == "breast") manifest <- manifest[manifest$cancer == "breast", , drop = FALSE]
  for (ii in seq_len(nrow(manifest))) {
    sample_id <- manifest$sample[[ii]]
    z <- res[res$cohort == cohort & res$sample == sample_id, , drop = FALSE]
    if (!nrow(z)) next
    z$rank_score <- ave(z$p_value, z$pathway, FUN = min)
    q_by_path <- tapply(z$q_value_sample_stat, z$pathway, min, na.rm = TRUE)
    p_by_path <- tapply(z$p_value, z$pathway, min, na.rm = TRUE)
    sig <- names(q_by_path)[q_by_path <= 0.10]
    top <- names(sort(p_by_path, decreasing = FALSE))[seq_len(min(max_pathways, length(p_by_path)))]
    chosen <- unique(c(sig, top))
    chosen <- chosen[order(p_by_path[chosen], q_by_path[chosen])]
    chosen <- head(chosen, max_pathways)
    sample_rows[[length(sample_rows) + 1L]] <- data.frame(
      cohort = cohort,
      cohort_task_id = ii,
      sample = sample_id,
      file = manifest$file[[ii]],
      selected_pathways = paste(chosen, collapse = ";"),
      n_selected_pathways = length(chosen),
      stringsAsFactors = FALSE
    )
  }
}

out <- do.call(rbind, sample_rows)
out$driver_task_id <- seq_len(nrow(out))
out <- out[, c("driver_task_id", "cohort", "cohort_task_id", "sample", "file", "n_selected_pathways", "selected_pathways")]
dir.create(file.path(root, "metadata"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(root, "metadata", "driver_stability_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# PASSAGE Driver Stability Manifest",
  "",
  paste0("- samples: ", nrow(out)),
  paste0("- max pathways per sample: ", max_pathways),
  paste0("- source pathway root: ", pathway_root)
), file.path(root, "metadata", "driver_stability_manifest.md"))
message("Wrote ", file.path(root, "metadata", "driver_stability_manifest.csv"))
