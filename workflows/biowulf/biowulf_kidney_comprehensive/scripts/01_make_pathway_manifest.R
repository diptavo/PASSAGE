#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
chunk_size <- if (length(args) >= 2L) as.integer(args[[2L]]) else 300L
kidney_root <- if (length(args) >= 3L) args[[3L]] else "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803"

pathway_file <- file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds")
meta_file <- file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered_metadata.csv")
if (!file.exists(pathway_file)) stop("Missing ", pathway_file)
pathways <- readRDS(pathway_file)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
manifest <- read.csv(file.path(kidney_root, "data", "passage_inputs", "passage_input_manifest.csv"), stringsAsFactors = FALSE)

tasks <- list()
tt <- 0L
for (ii in seq_len(nrow(manifest))) {
  for (coll in sort(unique(meta$collection))) {
    pws <- meta$pathway[meta$collection == coll]
    pws <- pws[pws %in% names(pathways)]
    if (!length(pws)) next
    chunks <- split(pws, ceiling(seq_along(pws) / chunk_size))
    for (cc in seq_along(chunks)) {
      tt <- tt + 1L
      tasks[[tt]] <- data.frame(
        task_id = tt,
        sample_task_id = ii,
        sample = manifest$sample[[ii]],
        file = manifest$file[[ii]],
        collection = coll,
        chunk = cc,
        n_pathways = length(chunks[[cc]]),
        pathways = paste(chunks[[cc]], collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}
out_dir <- file.path(root, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- do.call(rbind, tasks)
write.csv(out, file.path(out_dir, "pathway_task_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# Kidney MSigDB PASSAGE Pathway Manifest",
  "",
  paste0("- kidney samples/references: ", nrow(manifest)),
  paste0("- collections: ", paste(sort(unique(out$collection)), collapse = ", ")),
  paste0("- pathway tasks: ", nrow(out)),
  paste0("- chunk_size: ", chunk_size)
), file.path(out_dir, "pathway_task_manifest.md"))
message("Wrote ", nrow(out), " pathway tasks")
