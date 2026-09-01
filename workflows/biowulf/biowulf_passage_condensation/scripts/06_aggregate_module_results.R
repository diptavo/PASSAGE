#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(root, "results", "module_testing")
out_dir <- file.path(root, "results", "module_testing_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(res_root, pattern = "^module_task_.*[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No module result files under ", res_root)
tbl <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
tbl$q_value_sample_stat <- ave(tbl$p_value, tbl$sample, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
write.csv(tbl, file.path(out_dir, "module_all_results.csv"), row.names = FALSE)

fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}
meta <- do.call(rbind, lapply(split(tbl, interaction(tbl$cohort, tbl$reference, tbl$module_id, tbl$statistic, drop = TRUE, sep = "\r")), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    reference = z$reference[[1L]],
    module_id = z$module_id[[1L]],
    source_statistic = z$source_statistic[[1L]],
    statistic = z$statistic[[1L]],
    representative_pathway = z$representative_pathway[[1L]],
    n_member_pathways = z$n_member_pathways[[1L]],
    module_size_median = stats::median(z$module_size, na.rm = TRUE),
    n_samples = length(unique(z$sample)),
    fisher_p = fisher_p(z$p_value),
    min_p = min(z$p_value, na.rm = TRUE),
    median_effect_z = stats::median((z$observed - z$null_mean) / pmax(z$null_sd, .Machine$double.eps), na.rm = TRUE),
    frac_gpd_used = mean(z$gpd_used > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
meta$q_value_cohort_reference_stat <- ave(meta$fisher_p, meta$cohort, meta$reference, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta <- meta[order(meta$cohort, meta$reference, meta$statistic, meta$q_value_cohort_reference_stat, -meta$median_effect_z), , drop = FALSE]
write.csv(meta, file.path(out_dir, "module_reference_meta.csv"), row.names = FALSE)

diag <- read.csv(file.path(root, "results", "pathway_modules", "pathway_module_diagnostics.csv"), stringsAsFactors = FALSE)
top <- head(meta[order(meta$q_value_cohort_reference_stat, -meta$median_effect_z), , drop = FALSE], 80L)
md <- c(
  "# PASSAGE Module-Level Competitive Testing",
  "",
  paste0("- module result files: ", length(files)),
  paste0("- module result rows: ", nrow(tbl)),
  "",
  "## Compression Diagnostics",
  "",
  paste(c("cohort", "reference", "statistic", "significant_pathways", "modules", "compression"), collapse = " | "),
  paste(rep("---", 6), collapse = " | ")
)
for (ii in seq_len(nrow(diag))) {
  r <- diag[ii, ]
  md <- c(md, paste(c(r$cohort, r$reference, r$statistic, r$n_significant, r$n_modules, sprintf("%.2f", r$compression_ratio)), collapse = " | "))
}
md <- c(md, "", "## Top Module Signals", "",
        paste(c("cohort", "reference", "statistic", "module_id", "representative", "fisher_p", "q", "effect_z", "n_member_pathways"), collapse = " | "),
        paste(rep("---", 9), collapse = " | "))
for (ii in seq_len(nrow(top))) {
  r <- top[ii, ]
  md <- c(md, paste(c(r$cohort, r$reference, r$statistic, r$module_id, r$representative_pathway,
                      sprintf("%.4g", r$fisher_p), sprintf("%.4g", r$q_value_cohort_reference_stat),
                      sprintf("%.3f", r$median_effect_z), r$n_member_pathways), collapse = " | "))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
