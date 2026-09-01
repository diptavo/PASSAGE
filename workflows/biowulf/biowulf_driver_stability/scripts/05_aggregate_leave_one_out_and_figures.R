#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
fig_dir <- file.path(root, "results", "figures")
sum_dir <- file.path(root, "results", "driver_validation_summary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

loo_files <- list.files(file.path(root, "results", "leave_one_out"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
if (!length(loo_files)) stop("No leave-one-out CSV files found")
loo <- do.call(rbind, lapply(loo_files, read.csv, stringsAsFactors = FALSE))
write.csv(loo, file.path(sum_dir, "leave_one_out_all_results.csv"), row.names = FALSE)

drv_sum <- read.csv(file.path(root, "results", "driver_stability_summary", "driver_stability_by_cohort_statistic.csv"), stringsAsFactors = FALSE)
pair_sum <- read.csv(file.path(root, "results", "driver_stability_summary", "driver_stat_pairwise_topk_overlap_summary.csv"), stringsAsFactors = FALSE)
gene_freq <- read.csv(file.path(root, "results", "driver_stability_summary", "driver_stability_gene_frequencies.csv"), stringsAsFactors = FALSE)

loo_sum <- do.call(rbind, lapply(split(loo, list(loo$cohort, loo$statistic, loo$type), drop = TRUE), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    statistic = z$statistic[[1L]],
    type = z$type[[1L]],
    n = nrow(z),
    median_drop = stats::median(z$stat_drop, na.rm = TRUE),
    mean_drop = mean(z$stat_drop, na.rm = TRUE),
    frac_positive_drop = mean(z$stat_drop > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(loo_sum, file.path(sum_dir, "leave_one_out_by_cohort_statistic_type.csv"), row.names = FALSE)

driver_vs_control <- do.call(rbind, lapply(split(loo, list(loo$cohort, loo$statistic), drop = TRUE), function(z) {
  d <- z$stat_drop[z$type == "driver"]
  c <- z$stat_drop[z$type == "control"]
  wt <- suppressWarnings(stats::wilcox.test(d, c, alternative = "greater"))
  data.frame(
    cohort = z$cohort[[1L]],
    statistic = z$statistic[[1L]],
    n_driver = length(d),
    n_control = length(c),
    median_driver_drop = stats::median(d, na.rm = TRUE),
    median_control_drop = stats::median(c, na.rm = TRUE),
    median_delta = stats::median(d, na.rm = TRUE) - stats::median(c, na.rm = TRUE),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  )
}))
driver_vs_control$q_value <- ave(driver_vs_control$wilcox_p, driver_vs_control$cohort, FUN = function(x) p.adjust(x, "BH"))
write.csv(driver_vs_control, file.path(sum_dir, "leave_one_out_driver_vs_control.csv"), row.names = FALSE)

png(file.path(fig_dir, "01_driver_stability_jaccard_by_stat.png"), width = 1500, height = 900, res = 130)
op <- par(mar = c(9, 5, 3, 1))
labs <- paste(drv_sum$cohort, drv_sum$statistic, sep = "\n")
barplot(drv_sum$median_topk_jaccard, names.arg = labs, las = 2, col = "steelblue",
        ylab = "Median top-k Jaccard", main = "Driver Stability Across Statistics")
par(op)
dev.off()

png(file.path(fig_dir, "02_leave_one_out_driver_vs_control.png"), width = 1500, height = 900, res = 130)
op <- par(mar = c(9, 5, 3, 1))
plot_df <- driver_vs_control
labs <- paste(plot_df$cohort, plot_df$statistic, sep = "\n")
barplot(plot_df$median_delta, names.arg = labs, las = 2, col = ifelse(plot_df$median_delta > 0, "darkgreen", "firebrick"),
        ylab = "Median driver drop - median control drop", main = "Leave-One-Gene-Out Driver Validation")
abline(h = 0, lty = 2)
par(op)
dev.off()

png(file.path(fig_dir, "03_pairwise_driver_jaccard_vs_scorez.png"), width = 1400, height = 850, res = 130)
op <- par(mar = c(9, 5, 3, 1))
z <- pair_sum[pair_sum$stat2 %in% c("score_z", "score_z_robust") & !(pair_sum$stat1 %in% c("score_z", "score_z_robust")), , drop = FALSE]
labs <- paste(z$cohort, z$stat1, "vs", z$stat2, sep = "\n")
barplot(z$median_jaccard, names.arg = labs, las = 2, col = "darkorange",
        ylab = "Median top-k Jaccard", main = "Agreement With score_z Baselines")
par(op)
dev.off()

top_ex <- gene_freq[order(-gene_freq$selection_frequency, gene_freq$mean_rank), , drop = FALSE]
top_ex <- top_ex[top_ex$cohort %in% c("breast", "kidney", "dlpfc") & top_ex$statistic %in% c("CSPS", "GSPS", "score_z"), , drop = FALSE]
top_ex$key <- paste(top_ex$cohort, top_ex$sample, top_ex$pathway, top_ex$statistic, sep = "\r")
keys <- unique(top_ex$key)[seq_len(min(6L, length(unique(top_ex$key))))]
png(file.path(fig_dir, "04_representative_driver_frequency_bars.png"), width = 1500, height = 1000, res = 130)
op <- par(mfrow = c(2, 3), mar = c(7, 4, 3, 1))
for (kk in keys) {
  z <- head(top_ex[top_ex$key == kk, , drop = FALSE], 10L)
  barplot(z$selection_frequency, names.arg = z$gene, las = 2, col = "slateblue",
          ylim = c(0, 1), ylab = "Selection frequency",
          main = paste(z$cohort[[1L]], z$statistic[[1L]], z$pathway[[1L]], sep = "\n"))
}
par(op)
dev.off()

md <- c(
  "# PASSAGE Driver Validation",
  "",
  paste0("- leave-one-out rows: ", nrow(loo)),
  paste0("- figures: ", fig_dir),
  "",
  "## Leave-One-Out Driver vs Control",
  "",
  paste(c("cohort", "statistic", "median_driver_drop", "median_control_drop", "median_delta", "q"), collapse = " | "),
  paste(rep("---", 6), collapse = " | ")
)
driver_vs_control <- driver_vs_control[order(driver_vs_control$cohort, driver_vs_control$statistic), , drop = FALSE]
for (ii in seq_len(nrow(driver_vs_control))) {
  r <- driver_vs_control[ii, ]
  md <- c(md, paste(c(r$cohort, r$statistic,
                      sprintf("%.4g", r$median_driver_drop),
                      sprintf("%.4g", r$median_control_drop),
                      sprintf("%.4g", r$median_delta),
                      sprintf("%.4g", r$q_value)), collapse = " | "))
}
writeLines(md, file.path(sum_dir, "summary.md"))
message("Wrote driver validation outputs to ", sum_dir)
