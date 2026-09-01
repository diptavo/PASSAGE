#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
bench_root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(bench_root, "results", "pathway_new_stats")
out_dir <- file.path(bench_root, "results", "pathway_new_stats_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(res_root, pattern = "^pathway_new_stats_.*[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No result CSV files found under ", res_root)
tbl <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
tbl$q_value_sample_stat <- ave(tbl$p_value, tbl$sample, tbl$statistic, FUN = function(x) p.adjust(x, "BH"))
write.csv(tbl, file.path(out_dir, "pathway_new_stats_all_results.csv"), row.names = FALSE)

fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

split_key <- interaction(tbl$cohort, tbl$statistic, tbl$pathway, drop = TRUE, sep = "\r")
meta <- do.call(rbind, lapply(split(tbl, split_key), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    statistic = z$statistic[[1L]],
    pathway = z$pathway[[1L]],
    pathway_size_median = stats::median(z$pathway_size, na.rm = TRUE),
    n_samples = length(unique(z$sample)),
    fisher_p = fisher_p(z$p_value),
    min_p = min(z$p_value, na.rm = TRUE),
    median_p = stats::median(z$p_value, na.rm = TRUE),
    n_sample_q05 = sum(z$q_value_sample_stat <= 0.05, na.rm = TRUE),
    n_sample_q10 = sum(z$q_value_sample_stat <= 0.10, na.rm = TRUE),
    median_observed = stats::median(z$observed, na.rm = TRUE),
    median_null_mean = stats::median(z$null_mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
meta$q_value_cohort_stat <- ave(meta$fisher_p, meta$cohort, meta$statistic, FUN = function(x) p.adjust(x, "BH"))
meta <- meta[order(meta$cohort, meta$statistic, meta$q_value_cohort_stat, meta$fisher_p), , drop = FALSE]
write.csv(meta, file.path(out_dir, "pathway_new_stats_cohort_meta.csv"), row.names = FALSE)

sample_summary <- do.call(rbind, lapply(split(tbl, list(tbl$cohort, tbl$sample, tbl$statistic), drop = TRUE), function(z) {
  z <- z[order(z$q_value_sample_stat, z$p_value), , drop = FALSE]
  data.frame(
    cohort = z$cohort[[1L]],
    sample = z$sample[[1L]],
    statistic = z$statistic[[1L]],
    n_pathways = length(unique(z$pathway)),
    n_q05 = sum(z$q_value_sample_stat <= 0.05, na.rm = TRUE),
    n_q10 = sum(z$q_value_sample_stat <= 0.10, na.rm = TRUE),
    top_pathway = z$pathway[[1L]],
    top_p = z$p_value[[1L]],
    top_q = z$q_value_sample_stat[[1L]],
    stringsAsFactors = FALSE
  )
}))
sample_summary <- sample_summary[order(sample_summary$cohort, sample_summary$sample, sample_summary$statistic), , drop = FALSE]
write.csv(sample_summary, file.path(out_dir, "pathway_new_stats_sample_summary.csv"), row.names = FALSE)

top_meta <- meta[order(meta$q_value_cohort_stat, meta$fisher_p), , drop = FALSE]
top_meta <- head(top_meta, 40L)

md <- c(
  "# PASSAGE Pathway Analysis With New Spatial-Program Statistics",
  "",
  paste0("- result files: ", length(files)),
  paste0("- rows: ", nrow(tbl)),
  paste0("- cohorts: ", paste(sort(unique(tbl$cohort)), collapse = ", ")),
  paste0("- statistics: ", paste(sort(unique(tbl$statistic)), collapse = ", ")),
  paste0("- B: ", paste(sort(unique(tbl$B)), collapse = ", ")),
  "",
  "## Top Cohort-Level Pathways",
  "",
  paste(c("cohort", "statistic", "pathway", "fisher_p", "q", "n_samples", "n_sample_q05"), collapse = " | "),
  paste(rep("---", 7), collapse = " | ")
)
if (nrow(top_meta)) {
  for (ii in seq_len(nrow(top_meta))) {
    r <- top_meta[ii, ]
    md <- c(md, paste(c(r$cohort, r$statistic, r$pathway,
                        sprintf("%.4g", r$fisher_p), sprintf("%.4g", r$q_value_cohort_stat),
                        r$n_samples, r$n_sample_q05), collapse = " | "))
  }
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
