#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(root, "results", "pathway_testing")
out_dir <- file.path(root, "results", "pathway_testing_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(res_root, pattern = "^pathway_task_.*[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No pathway result files under ", res_root)
tbl <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
tbl$q_value_sample_stat <- ave(tbl$p_value, tbl$sample, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
tbl$q_value_sample_collection_stat <- ave(tbl$p_value, tbl$sample, tbl$collection, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
write.csv(tbl, file.path(out_dir, "kidney_msigdb_pathway_all_results.csv"), row.names = FALSE)

fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

split_key <- interaction(tbl$reference, tbl$collection, tbl$statistic, tbl$pathway, drop = TRUE, sep = "\r")
meta <- do.call(rbind, lapply(split(tbl, split_key), function(z) {
  data.frame(
    reference = z$reference[[1L]],
    collection = z$collection[[1L]],
    statistic = z$statistic[[1L]],
    pathway = z$pathway[[1L]],
    pathway_size_median = stats::median(z$pathway_size, na.rm = TRUE),
    n_samples = length(unique(z$sample)),
    fisher_p = fisher_p(z$p_value),
    min_p = min(z$p_value, na.rm = TRUE),
    median_p = stats::median(z$p_value, na.rm = TRUE),
    n_sample_q05 = sum(z$q_value_sample_stat <= 0.05, na.rm = TRUE),
    n_sample_collection_q05 = sum(z$q_value_sample_collection_stat <= 0.05, na.rm = TRUE),
    median_observed = stats::median(z$observed, na.rm = TRUE),
    median_null_mean = stats::median(z$null_mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
meta$q_value_reference_stat <- ave(meta$fisher_p, meta$reference, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta$q_value_reference_collection_stat <- ave(meta$fisher_p, meta$reference, meta$collection, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta <- meta[order(meta$reference, meta$statistic, meta$q_value_reference_stat, meta$fisher_p), , drop = FALSE]
write.csv(meta, file.path(out_dir, "kidney_msigdb_pathway_reference_meta.csv"), row.names = FALSE)

sample_summary <- do.call(rbind, lapply(split(tbl, list(tbl$sample, tbl$collection, tbl$statistic), drop = TRUE), function(z) {
  z <- z[order(z$q_value_sample_collection_stat, z$p_value), , drop = FALSE]
  data.frame(
    sample = z$sample[[1L]],
    reference = z$reference[[1L]],
    collection = z$collection[[1L]],
    statistic = z$statistic[[1L]],
    n_pathways = length(unique(z$pathway)),
    n_q05 = sum(z$q_value_sample_collection_stat <= 0.05, na.rm = TRUE),
    n_q10 = sum(z$q_value_sample_collection_stat <= 0.10, na.rm = TRUE),
    top_pathway = z$pathway[[1L]],
    top_p = z$p_value[[1L]],
    top_q = z$q_value_sample_collection_stat[[1L]],
    stringsAsFactors = FALSE
  )
}))
write.csv(sample_summary, file.path(out_dir, "kidney_msigdb_pathway_sample_summary.csv"), row.names = FALSE)

top_meta <- head(meta[order(meta$q_value_reference_stat, meta$fisher_p), , drop = FALSE], 80L)
md <- c(
  "# Comprehensive Kidney PASSAGE MSigDB Pathway Testing",
  "",
  paste0("- result_files: ", length(files)),
  paste0("- rows: ", nrow(tbl)),
  paste0("- samples/references: ", length(unique(tbl$sample))),
  paste0("- collections: ", paste(sort(unique(tbl$collection)), collapse = ", ")),
  paste0("- statistics: ", paste(sort(unique(tbl$statistic)), collapse = ", ")),
  "",
  "## Top Reference-Level Signals",
  "",
  paste(c("reference", "collection", "statistic", "pathway", "fisher_p", "q_reference_stat", "n_samples", "n_sample_q05"), collapse = " | "),
  paste(rep("---", 8), collapse = " | ")
)
for (ii in seq_len(nrow(top_meta))) {
  r <- top_meta[ii, ]
  md <- c(md, paste(c(r$reference, r$collection, r$statistic, r$pathway,
                      sprintf("%.4g", r$fisher_p), sprintf("%.4g", r$q_value_reference_stat),
                      r$n_samples, r$n_sample_q05), collapse = " | "))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
