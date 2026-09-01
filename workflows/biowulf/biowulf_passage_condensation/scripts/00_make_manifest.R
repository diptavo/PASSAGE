#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
cohorts <- if (length(args) >= 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("kidney", "breast", "dlpfc")
chunk_size <- if (length(args) >= 3L) as.integer(args[[3L]]) else 300L

data_roots <- c(
  kidney = "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803",
  breast = "/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803",
  dlpfc = "/data/Dutta_lab/SPATH/PASSAGE_spatialDLPFC_20260802"
)

pathway_file <- file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds")
meta_file <- file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered_metadata.csv")
if (!file.exists(pathway_file)) stop("Missing ", pathway_file)
pathways <- readRDS(pathway_file)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

read_cohort_manifest <- function(cohort) {
  m <- read.csv(file.path(data_roots[[cohort]], "data", "passage_inputs", "passage_input_manifest.csv"), stringsAsFactors = FALSE)
  if (cohort == "breast" && "cancer" %in% names(m)) m <- m[m$cancer == "breast", , drop = FALSE]
  if (!"reference" %in% names(m)) m$reference <- paste0(cohort, "_celltype_broad")
  if (!"spatial_sample" %in% names(m)) m$spatial_sample <- m$sample
  m$cohort <- cohort
  m[, c("cohort", "sample", "spatial_sample", "reference", "file"), drop = FALSE]
}

sample_tbl <- do.call(rbind, lapply(cohorts, read_cohort_manifest))
tasks <- list()
tt <- 0L
for (ii in seq_len(nrow(sample_tbl))) {
  for (coll in sort(unique(meta$collection))) {
    pws <- meta$pathway[meta$collection == coll]
    pws <- pws[pws %in% names(pathways)]
    chunks <- split(pws, ceiling(seq_along(pws) / chunk_size))
    for (cc in seq_along(chunks)) {
      tt <- tt + 1L
      tasks[[tt]] <- data.frame(
        task_id = tt,
        cohort = sample_tbl$cohort[[ii]],
        sample_task_id = ii,
        sample = sample_tbl$sample[[ii]],
        spatial_sample = sample_tbl$spatial_sample[[ii]],
        reference = sample_tbl$reference[[ii]],
        file = sample_tbl$file[[ii]],
        collection = coll,
        chunk = cc,
        n_pathways = length(chunks[[cc]]),
        pathways = paste(chunks[[cc]], collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}
out <- do.call(rbind, tasks)
dir.create(file.path(root, "metadata"), recursive = TRUE, showWarnings = FALSE)
write.csv(sample_tbl, file.path(root, "metadata", "sample_manifest.csv"), row.names = FALSE)
write.csv(out, file.path(root, "metadata", "pathway_task_manifest.csv"), row.names = FALSE)
writeLines(c(
  "# PASSAGE Condensation Manifest",
  "",
  paste0("- cohorts: ", paste(cohorts, collapse = ", ")),
  paste0("- samples/references: ", nrow(sample_tbl)),
  paste0("- pathway tasks: ", nrow(out)),
  paste0("- chunk_size: ", chunk_size),
  paste0("- collections: ", paste(sort(unique(out$collection)), collapse = ", "))
), file.path(root, "metadata", "manifest_summary.md"))
message("Wrote ", nrow(out), " pathway tasks")
