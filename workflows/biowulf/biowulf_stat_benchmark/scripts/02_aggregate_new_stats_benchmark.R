#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
bench_root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(bench_root, "results", "new_stats_benchmark")
out_dir <- file.path(bench_root, "results", "new_stats_benchmark_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(res_root, pattern = "^new_stats_.*[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No result CSV files found under ", res_root)
tbl <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
write.csv(tbl, file.path(out_dir, "new_stats_all_pvalues.csv"), row.names = FALSE)

alpha_grid <- c(0.10, 0.05, 0.01)
summary <- do.call(rbind, lapply(split(tbl, list(tbl$cohort, tbl$statistic), drop = TRUE), function(z) {
  do.call(rbind, lapply(alpha_grid, function(a) {
    data.frame(
      cohort = z$cohort[[1L]],
      statistic = z$statistic[[1L]],
      alpha = a,
      n = nrow(z),
      n_samples = length(unique(z$sample)),
      reject_rate = mean(z$p_value <= a, na.rm = TRUE),
      mc_se = sqrt(pmax(mean(z$p_value <= a, na.rm = TRUE) * (1 - mean(z$p_value <= a, na.rm = TRUE)), 0) / nrow(z)),
      median_p = stats::median(z$p_value, na.rm = TRUE),
      min_p = min(z$p_value, na.rm = TRUE),
      median_elapsed_sec = stats::median(z$elapsed_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
summary <- summary[order(summary$cohort, summary$statistic, summary$alpha), , drop = FALSE]
write.csv(summary, file.path(out_dir, "new_stats_type1_summary.csv"), row.names = FALSE)

overall <- do.call(rbind, lapply(split(tbl, tbl$statistic), function(z) {
  do.call(rbind, lapply(alpha_grid, function(a) {
    data.frame(
      cohort = "ALL",
      statistic = z$statistic[[1L]],
      alpha = a,
      n = nrow(z),
      n_samples = length(unique(paste(z$cohort, z$sample, sep = "::"))),
      reject_rate = mean(z$p_value <= a, na.rm = TRUE),
      mc_se = sqrt(pmax(mean(z$p_value <= a, na.rm = TRUE) * (1 - mean(z$p_value <= a, na.rm = TRUE)), 0) / nrow(z)),
      median_p = stats::median(z$p_value, na.rm = TRUE),
      min_p = min(z$p_value, na.rm = TRUE),
      median_elapsed_sec = stats::median(z$elapsed_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
overall <- overall[order(overall$statistic, overall$alpha), , drop = FALSE]
write.csv(overall, file.path(out_dir, "new_stats_type1_summary_overall.csv"), row.names = FALSE)

md <- c(
  "# PASSAGE New Spatial-Program Statistic Benchmark",
  "",
  paste0("- result files: ", length(files)),
  paste0("- rows: ", nrow(tbl)),
  paste0("- cohorts: ", paste(sort(unique(tbl$cohort)), collapse = ", ")),
  "",
  "## Overall Type I Summary",
  "",
  paste(c("statistic", "alpha", "n", "n_samples", "reject_rate", "mc_se", "median_p", "min_p"), collapse = " | "),
  paste(rep("---", 8), collapse = " | ")
)
for (ii in seq_len(nrow(overall))) {
  r <- overall[ii, ]
  md <- c(md, paste(c(r$statistic, r$alpha, r$n, r$n_samples,
                      sprintf("%.4f", r$reject_rate), sprintf("%.4f", r$mc_se),
                      sprintf("%.4f", r$median_p), sprintf("%.4f", r$min_p)),
                    collapse = " | "))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
