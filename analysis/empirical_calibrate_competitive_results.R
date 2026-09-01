# Empirically calibrate saved PASSAGE competitive residual-null p-values.
#
# Usage:
#   Rscript scripts/empirical_calibrate_competitive_results.R \
#     --result-dirs=results/run_a,results/run_b \
#     --out-dir=results/empirical_competitive_calibration

parse_args <- function(args) {
  cfg <- list(
    result_dirs = character(),
    out_dir = file.path("results", "empirical_competitive_calibration")
  )
  for (arg in args) {
    if (grepl("^--result-dirs=", arg)) {
      cfg$result_dirs <- strsplit(sub("^--result-dirs=", "", arg), ",", fixed = TRUE)[[1]]
    }
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
  }
  if (!length(cfg$result_dirs)) stop("--result-dirs is required")
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

summarize_endpoint <- function(x, p_col) {
  groups <- unique(x[c("dataset", "mode", "matching", "mc_sampler")])
  rows <- list()
  ii <- 0L
  for (gg in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(x))
    for (cc in colnames(groups)) keep <- keep & x[[cc]] == groups[[cc]][gg]
    z <- x[keep, , drop = FALSE]
    for (alpha in c(0.10, 0.05, 0.01)) {
      p <- z[[p_col]][is.finite(z[[p_col]])]
      rate <- mean(p <= alpha)
      n <- length(p)
      se <- sqrt(rate * (1 - rate) / max(1L, n))
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        groups[gg, , drop = FALSE],
        endpoint = p_col,
        alpha = alpha,
        n = n,
        reject_rate = rate,
        mc_se = se,
        ci95_low = max(0, rate - 1.96 * se),
        ci95_high = min(1, rate + 1.96 * se),
        median_p = if (length(p)) stats::median(p) else NA_real_,
        min_p = if (length(p)) min(p) else NA_real_,
        median_elapsed_sec = stats::median(z$elapsed_sec),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

pieces <- lapply(cfg$result_dirs, function(d) {
  f <- file.path(d, "real_residual_null_competitive_pvalues.csv")
  if (!file.exists(f)) stop("missing p-value file: ", f)
  x <- read.csv(f, stringsAsFactors = FALSE)
  x$calibration_run <- basename(normalizePath(d, mustWork = FALSE))
  x
})
out <- do.call(rbind, pieces)
out <- passage_empirical_competitive_calibration(out)

for (p_col in c("competitive_score_p",
                "competitive_score_empirical_pooled_p",
                "competitive_score_empirical_leave_rep_p",
                "competitive_score_empirical_size_p",
                "competitive_score_empirical_leave_rep_size_p")) {
  fdr_col <- sub("_p$", "_fdr", p_col)
  out[[fdr_col]] <- stats::p.adjust(out[[p_col]], method = "BH")
}

summary <- rbind(
  summarize_endpoint(out, "competitive_score_p"),
  summarize_endpoint(out, "competitive_score_empirical_pooled_p"),
  summarize_endpoint(out, "competitive_score_empirical_leave_rep_p"),
  summarize_endpoint(out, "competitive_score_empirical_size_p"),
  summarize_endpoint(out, "competitive_score_empirical_leave_rep_size_p")
)

write.csv(out, file.path(cfg$out_dir, "empirical_calibrated_competitive_pvalues.csv"), row.names = FALSE)
write.csv(summary, file.path(cfg$out_dir, "empirical_calibrated_competitive_summary.csv"), row.names = FALSE)

md <- c(
  "# PASSAGE Competitive Empirical Calibration",
  "",
  paste0("- Result directories: ", paste(cfg$result_dirs, collapse = ", ")),
  "",
  paste(colnames(summary), collapse = " | "),
  paste(rep("---", ncol(summary)), collapse = " | "),
  apply(summary, 1L, function(x) paste(x, collapse = " | "))
)
writeLines(md, file.path(cfg$out_dir, "summary.md"))
message("Wrote outputs to ", cfg$out_dir)
print(summary, row.names = FALSE)
