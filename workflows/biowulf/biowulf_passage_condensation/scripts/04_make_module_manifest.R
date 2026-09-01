#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
max_modules_per_group <- if (length(args) >= 2L) as.integer(args[[2L]]) else 250L

samples <- read.csv(file.path(root, "metadata", "sample_manifest.csv"), stringsAsFactors = FALSE)
modules <- read.csv(file.path(root, "results", "pathway_modules", "pathway_modules.csv"), stringsAsFactors = FALSE)
if (!nrow(modules)) stop("No modules found")
modules <- modules[order(modules$cohort, modules$reference, modules$source_statistic, modules$best_q, -modules$median_effect_z), , drop = FALSE]
modules <- do.call(rbind, lapply(split(modules, list(modules$cohort, modules$reference, modules$source_statistic), drop = TRUE), head, n = max_modules_per_group))

tasks <- list()
tt <- 0L
for (ii in seq_len(nrow(samples))) {
  z <- modules[modules$cohort == samples$cohort[[ii]] & modules$reference == samples$reference[[ii]], , drop = FALSE]
  if (!nrow(z)) next
  for (jj in seq_len(nrow(z))) {
    tt <- tt + 1L
    tasks[[tt]] <- data.frame(
      task_id = tt,
      cohort = samples$cohort[[ii]],
      sample = samples$sample[[ii]],
      spatial_sample = samples$spatial_sample[[ii]],
      reference = samples$reference[[ii]],
      file = samples$file[[ii]],
      module_id = z$module_id[[jj]],
      source_statistic = z$source_statistic[[jj]],
      representative_pathway = z$representative_pathway[[jj]],
      n_member_pathways = z$n_member_pathways[[jj]],
      n_module_genes = z$n_module_genes[[jj]],
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, tasks)
write.csv(out, file.path(root, "metadata", "module_task_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# PASSAGE Module Test Manifest",
  "",
  paste0("- module tasks: ", nrow(out)),
  paste0("- max_modules_per_group: ", max_modules_per_group)
), file.path(root, "metadata", "module_task_manifest.md"))
message("Wrote ", nrow(out), " module tasks")
