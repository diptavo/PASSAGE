#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()

cache_dir <- file.path(root, "refs", "msigdb", "r_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = cache_dir)
options(timeout = max(600, getOption("timeout")))

ensure_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  install.packages(pkg, repos = "https://cloud.r-project.org")
  invisible(TRUE)
}

ensure_pkg("msigdbr")

out_dir <- file.path(root, "refs", "msigdb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "hallmark_human_pathways.rds")

hall <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hall$gene_symbol <- toupper(hall$gene_symbol)
pathways <- lapply(split(hall$gene_symbol, hall$gs_name), unique)
saveRDS(pathways, out_file)

summary_file <- file.path(out_dir, "hallmark_human_pathways_summary.txt")
writeLines(c(
  paste0("file: ", out_file),
  paste0("created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("n_pathways: ", length(pathways)),
  paste0("n_unique_genes: ", length(unique(unlist(pathways, use.names = FALSE)))),
  paste0("msigdbr_version: ", as.character(utils::packageVersion("msigdbr")))
), summary_file)

message("Wrote ", out_file)
