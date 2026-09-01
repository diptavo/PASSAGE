# Compare two 10x PASSAGE Hallmark result roots, typically the primary
# technical-covariate run and a spatial-covariate sensitivity run.
#
# Usage:
#   Rscript scripts/compare_10x_passage_sensitivity.R \
#     --baseline-root=results/passage_10x_hallmark_perm999_fastcal \
#     --sensitivity-root=results/passage_10x_hallmark_spatial_quad_perm999_fastcal

parse_args <- function(args) {
  cfg <- list(
    baseline_root = "results/passage_10x_hallmark_perm999_fastcal",
    sensitivity_root = "results/passage_10x_hallmark_spatial_quad_perm999_fastcal",
    out_file = NULL
  )
  for (arg in args) {
    if (grepl("^--baseline-root=", arg)) cfg$baseline_root <- sub("^--baseline-root=", "", arg)
    if (grepl("^--sensitivity-root=", arg)) cfg$sensitivity_root <- sub("^--sensitivity-root=", "", arg)
    if (grepl("^--out-file=", arg)) cfg$out_file <- sub("^--out-file=", "", arg)
  }
  if (is.null(cfg$out_file)) {
    cfg$out_file <- file.path(cfg$sensitivity_root, "sensitivity_comparison.csv")
  }
  cfg
}

read_result <- function(root, dataset) {
  read.csv(file.path(root, dataset, "passage_hallmark_pathways.csv"), stringsAsFactors = FALSE)
}

compare_dataset <- function(dataset, cfg) {
  a <- read_result(cfg$baseline_root, dataset)
  b <- read_result(cfg$sensitivity_root, dataset)
  m <- merge(
    a[, c("pathway", "p_H1", "fdr_H1", "R2_cca", "PSVS_range", "spasset_size", "spasset_genes")],
    b[, c("pathway", "p_H1", "fdr_H1", "R2_cca", "PSVS_range", "spasset_size", "spasset_genes")],
    by = "pathway",
    suffixes = c("_baseline", "_sensitivity")
  )
  m$dataset <- dataset
  m$delta_R2 <- m$R2_cca_sensitivity - m$R2_cca_baseline
  m$delta_PSVS_range <- m$PSVS_range_sensitivity - m$PSVS_range_baseline
  m$significant_baseline <- m$fdr_H1_baseline <= 0.05
  m$significant_sensitivity <- m$fdr_H1_sensitivity <= 0.05
  m <- m[order(m$p_H1_sensitivity, -m$R2_cca_sensitivity, m$pathway), , drop = FALSE]
  rownames(m) <- NULL
  m
}

write_summary_md <- function(comparison, cfg) {
  out_md <- sub("[.]csv$", ".md", cfg$out_file)
  datasets <- unique(comparison$dataset)
  lines <- c(
    "# PASSAGE Hallmark Sensitivity Comparison",
    "",
    paste0("- Baseline root: `", cfg$baseline_root, "`"),
    paste0("- Sensitivity root: `", cfg$sensitivity_root, "`"),
    "",
    "dataset | baseline FDR<=0.05 | sensitivity FDR<=0.05 | median baseline R2 | median sensitivity R2",
    "--- | ---: | ---: | ---: | ---:"
  )
  for (dataset in datasets) {
    x <- comparison[comparison$dataset == dataset, , drop = FALSE]
    lines <- c(lines, paste(
      dataset,
      sum(x$significant_baseline),
      sum(x$significant_sensitivity),
      sprintf("%.3f", median(x$R2_cca_baseline)),
      sprintf("%.3f", median(x$R2_cca_sensitivity)),
      sep = " | "
    ))
  }
  lines <- c(lines, "", "## Top Sensitivity Hits", "")
  for (dataset in datasets) {
    x <- comparison[comparison$dataset == dataset, , drop = FALSE]
    top <- head(x[, c("pathway", "p_H1_sensitivity", "fdr_H1_sensitivity", "R2_cca_baseline", "R2_cca_sensitivity", "delta_R2")], 12)
    lines <- c(
      lines,
      paste0("### ", dataset),
      "",
      "pathway | p sensitivity | fdr sensitivity | R2 baseline | R2 sensitivity | delta R2",
      "--- | ---: | ---: | ---: | ---: | ---:"
    )
    rows <- apply(top, 1L, function(r) paste(r, collapse = " | "))
    lines <- c(lines, rows, "")
  }
  writeLines(lines, out_md)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
datasets <- intersect(
  list.dirs(cfg$baseline_root, recursive = FALSE, full.names = FALSE),
  list.dirs(cfg$sensitivity_root, recursive = FALSE, full.names = FALSE)
)
if (length(datasets) == 0L) {
  stop("No overlapping dataset folders found")
}
comparison <- do.call(rbind, lapply(datasets, compare_dataset, cfg = cfg))
dir.create(dirname(cfg$out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, cfg$out_file, row.names = FALSE)
write_summary_md(comparison, cfg)
message("Wrote ", cfg$out_file)
message("Wrote ", sub("[.]csv$", ".md", cfg$out_file))
