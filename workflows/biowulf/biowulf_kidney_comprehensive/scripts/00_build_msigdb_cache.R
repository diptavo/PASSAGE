#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
collections <- if (length(args) >= 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("H", "C2", "C3", "C5", "C6", "C7", "C8")
min_genes <- if (length(args) >= 3L) as.integer(args[[3L]]) else 15L
max_genes <- if (length(args) >= 4L) as.integer(args[[4L]]) else 250L

lib <- file.path(root, "rlib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr", repos = "https://cloud.r-project.org", lib = lib, quiet = TRUE)
}
if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Could not load or install msigdbr")

out_dir <- file.path(root, "refs", "msigdb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fetch_collection <- function(coll) {
  f <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% f) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = coll)
  } else if ("category" %in% f) {
    msigdbr::msigdbr(species = "Homo sapiens", category = coll)
  } else {
    x <- msigdbr::msigdbr(species = "Homo sapiens")
    x[x$gs_cat == coll | x$gs_collection == coll, , drop = FALSE]
  }
}

all_sets <- list()
meta_rows <- list()
for (coll in collections) {
  message("Fetching MSigDB collection ", coll)
  x <- fetch_collection(coll)
  if (!nrow(x)) next
  coll_col <- if ("gs_collection" %in% names(x)) "gs_collection" else "gs_cat"
  sub_col <- if ("gs_subcollection" %in% names(x)) "gs_subcollection" else if ("gs_subcat" %in% names(x)) "gs_subcat" else NA_character_
  split_sets <- split(toupper(x$gene_symbol), x$gs_name)
  for (nm in names(split_sets)) {
    genes <- sort(unique(split_sets[[nm]]))
    if (length(genes) < min_genes || length(genes) > max_genes) next
    key <- make.names(paste(coll, nm, sep = "__"), unique = TRUE)
    all_sets[[key]] <- genes
    one <- x[x$gs_name == nm, , drop = FALSE]
    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      pathway = key,
      msigdb_name = nm,
      collection = unique(as.character(one[[coll_col]]))[1L],
      subcollection = if (!is.na(sub_col)) unique(as.character(one[[sub_col]]))[1L] else "",
      n_genes = length(genes),
      stringsAsFactors = FALSE
    )
  }
}

if (!length(all_sets)) stop("No MSigDB sets passed filters")
meta <- do.call(rbind, meta_rows)
saveRDS(all_sets, file.path(out_dir, "msigdb_human_pathways_filtered.rds"))
write.csv(meta, file.path(out_dir, "msigdb_human_pathways_filtered_metadata.csv"), row.names = FALSE)
writeLines(c(
  "# MSigDB Cache",
  "",
  paste0("- collections: ", paste(collections, collapse = ", ")),
  paste0("- min_genes: ", min_genes),
  paste0("- max_genes: ", max_genes),
  paste0("- retained_sets: ", length(all_sets)),
  paste0("- msigdbr_version: ", as.character(utils::packageVersion("msigdbr")))
), file.path(out_dir, "README.md"))
message("Wrote ", length(all_sets), " filtered MSigDB gene sets")
