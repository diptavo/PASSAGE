#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(root, "results", "pathway_testing")
out_dir <- file.path(root, "results", "pathway_testing_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(res_root, pattern = "^pathway_task_.*[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No pathway result files under ", res_root)
tbl <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
tbl$p_value <- pmin(tbl$p_empirical, tbl$p_gpd, na.rm = TRUE)
tbl$q_value_sample_stat <- ave(tbl$p_value, tbl$sample, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
tbl$q_value_sample_collection_stat <- ave(tbl$p_value, tbl$sample, tbl$collection, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
write.csv(tbl, file.path(out_dir, "pathway_all_results.csv"), row.names = FALSE)

fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

split_key <- interaction(tbl$cohort, tbl$reference, tbl$collection, tbl$statistic, tbl$pathway, drop = TRUE, sep = "\r")
meta <- do.call(rbind, lapply(split(tbl, split_key), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    reference = z$reference[[1L]],
    collection = z$collection[[1L]],
    statistic = z$statistic[[1L]],
    pathway = z$pathway[[1L]],
    pathway_size_median = stats::median(z$pathway_size, na.rm = TRUE),
    n_samples = length(unique(z$sample)),
    fisher_p = fisher_p(z$p_value),
    fisher_p_empirical = fisher_p(z$p_empirical),
    min_p = min(z$p_value, na.rm = TRUE),
    median_p = stats::median(z$p_value, na.rm = TRUE),
    n_sample_q05 = sum(z$q_value_sample_stat <= 0.05, na.rm = TRUE),
    n_sample_collection_q05 = sum(z$q_value_sample_collection_stat <= 0.05, na.rm = TRUE),
    median_effect_z = stats::median((z$observed - z$null_mean) / pmax(z$null_sd, .Machine$double.eps), na.rm = TRUE),
    frac_gpd_used = mean(z$gpd_used > 0, na.rm = TRUE),
    genes = z$genes[[which.min(z$p_value)]],
    stringsAsFactors = FALSE
  )
}))
meta$q_value_cohort_reference_stat <- ave(meta$fisher_p, meta$cohort, meta$reference, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta$q_value_cohort_reference_collection_stat <- ave(meta$fisher_p, meta$cohort, meta$reference, meta$collection, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta <- meta[order(meta$cohort, meta$reference, meta$statistic, meta$q_value_cohort_reference_stat, meta$fisher_p), , drop = FALSE]
write.csv(meta, file.path(out_dir, "pathway_reference_meta.csv"), row.names = FALSE)

diag <- do.call(rbind, lapply(split(meta, list(meta$cohort, meta$reference, meta$statistic), drop = TRUE), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    reference = z$reference[[1L]],
    statistic = z$statistic[[1L]],
    n_tested_pathways = nrow(z),
    n_q05 = sum(z$q_value_cohort_reference_stat <= 0.05, na.rm = TRUE),
    n_q10 = sum(z$q_value_cohort_reference_stat <= 0.10, na.rm = TRUE),
    median_frac_gpd_used = stats::median(z$frac_gpd_used, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(diag, file.path(out_dir, "pathway_saturation_diagnostics.csv"), row.names = FALSE)

top <- head(meta[order(meta$q_value_cohort_reference_stat, -meta$median_effect_z), , drop = FALSE], 80L)
md <- c(
  "# PASSAGE Broad Pathway Results",
  "",
  paste0("- result files: ", length(files)),
  paste0("- rows: ", nrow(tbl)),
  paste0("- cohorts: ", paste(sort(unique(tbl$cohort)), collapse = ", ")),
  paste0("- statistics: ", paste(sort(unique(tbl$statistic)), collapse = ", ")),
  "",
  "## Saturation Diagnostics",
  "",
  paste(c("cohort", "reference", "statistic", "tested", "q05", "q10", "median_frac_gpd_used"), collapse = " | "),
  paste(rep("---", 7), collapse = " | ")
)
for (ii in seq_len(nrow(diag))) {
  r <- diag[ii, ]
  md <- c(md, paste(c(r$cohort, r$reference, r$statistic, r$n_tested_pathways, r$n_q05, r$n_q10, sprintf("%.3f", r$median_frac_gpd_used)), collapse = " | "))
}
md <- c(md, "", "## Top Signals", "",
        paste(c("cohort", "reference", "collection", "statistic", "pathway", "fisher_p", "q", "effect_z"), collapse = " | "),
        paste(rep("---", 8), collapse = " | "))
for (ii in seq_len(nrow(top))) {
  r <- top[ii, ]
  md <- c(md, paste(c(r$cohort, r$reference, r$collection, r$statistic, r$pathway,
                      sprintf("%.4g", r$fisher_p), sprintf("%.4g", r$q_value_cohort_reference_stat),
                      sprintf("%.3f", r$median_effect_z)), collapse = " | "))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
