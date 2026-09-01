#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)
source(file.path(root, "scripts", "load_passage.R"))

out_dir <- file.path(root, "results", "passage_cancer_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_many <- function(files) {
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

manifest <- read.csv(file.path(root, "data", "passage_inputs", "passage_input_manifest.csv"),
                     stringsAsFactors = FALSE)
h_files <- file.path(root, "results", "passage_cancer_h1_h3", manifest$sample,
                     "passage_cancer_h1_h3_pathways.csv")
c_files <- file.path(root, "results", "passage_cancer_competitive_h3", manifest$sample,
                     "passage_cancer_competitive_h3.csv")
h <- read_many(h_files)
c <- read_many(c_files)

combine_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}
split_key <- function(df, cols) interaction(df[, cols, drop = FALSE], drop = TRUE, sep = "||")

if (nrow(h)) {
  h_summary <- do.call(rbind, lapply(split(h, split_key(h, c("cancer", "reference", "pathway"))), function(z) {
    data.frame(
      cancer = z$cancer[[1L]],
      reference = z$reference[[1L]],
      pathway = z$pathway[[1L]],
      n_samples = length(unique(z$sample)),
      n_H1_fdr05 = sum(z$fdr_H1 <= 0.05, na.rm = TRUE),
      n_H3_fdr05 = sum(z$fdr_H3 <= 0.05, na.rm = TRUE),
      median_p_H1 = stats::median(z$p_H1, na.rm = TRUE),
      median_p_H3 = stats::median(z$p_H3, na.rm = TRUE),
      fisher_p_H1 = combine_p(z$p_H1),
      fisher_p_H3 = combine_p(z$p_H3),
      median_celltype_adjusted_reduction = stats::median(z$celltype_adjusted_reduction, na.rm = TRUE),
      median_H3_over_H1_Q = stats::median(z$H3_over_H1_Q, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  h_summary$fisher_fdr_H1 <- ave(h_summary$fisher_p_H1, h_summary$cancer, h_summary$reference,
                                 FUN = function(p) stats::p.adjust(p, method = "BH"))
  h_summary$fisher_fdr_H3 <- ave(h_summary$fisher_p_H3, h_summary$cancer, h_summary$reference,
                                 FUN = function(p) stats::p.adjust(p, method = "BH"))
  h_summary <- h_summary[order(h_summary$cancer, h_summary$reference,
                               h_summary$fisher_p_H3, h_summary$fisher_p_H1), , drop = FALSE]
} else {
  h_summary <- data.frame()
}

if (nrow(c)) {
  c$dataset <- paste(c$cancer, c$reference, sep = "_")
  c$mode <- "h3_marker_adjusted_score"
  c$matching <- "expr_detect_var_spatial"
  c$mc_sampler <- "module"
  c$replicate <- c$sample
  c <- passage_empirical_tail_calibrate(
    c, p_col = "competitive_score_p",
    group_cols = c("dataset", "statistic", "mode", "matching", "mc_sampler"),
    leaveout_cols = "spatial_sample",
    out_col = "competitive_score_empirical_leave_sample_p",
    n_col = "competitive_score_empirical_leave_sample_n"
  )
  c$competitive_score_fdr_global <- ave(c$competitive_score_p, c$cancer, c$reference, c$statistic,
                                        FUN = function(p) stats::p.adjust(p, method = "BH"))
  c$competitive_score_empirical_leave_sample_fdr <- ave(
    c$competitive_score_empirical_leave_sample_p, c$cancer, c$reference, c$statistic,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
  c_summary <- do.call(rbind, lapply(split(c, split_key(c, c("cancer", "reference", "statistic", "pathway"))), function(z) {
    data.frame(
      cancer = z$cancer[[1L]],
      reference = z$reference[[1L]],
      statistic = z$statistic[[1L]],
      pathway = z$pathway[[1L]],
      n_samples = length(unique(z$sample)),
      n_raw_p05 = sum(z$competitive_score_p <= 0.05, na.rm = TRUE),
      n_raw_fdr05 = sum(z$competitive_score_fdr_global <= 0.05, na.rm = TRUE),
      n_empirical_fdr05 = sum(z$competitive_score_empirical_leave_sample_fdr <= 0.05, na.rm = TRUE),
      median_p = stats::median(z$competitive_score_p, na.rm = TRUE),
      median_empirical_p = stats::median(z$competitive_score_empirical_leave_sample_p, na.rm = TRUE),
      fisher_p = combine_p(z$competitive_score_p),
      median_cEPSV = stats::median(z$cEPSV, na.rm = TRUE),
      median_coherence_pc1 = stats::median(z$coherence_pc1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  c_summary$fisher_fdr <- ave(c_summary$fisher_p, c_summary$cancer, c_summary$reference, c_summary$statistic,
                              FUN = function(p) stats::p.adjust(p, method = "BH"))
  c_summary <- c_summary[order(c_summary$cancer, c_summary$reference,
                               c_summary$statistic, c_summary$fisher_p), , drop = FALSE]
} else {
  c_summary <- data.frame()
}

write.csv(h, file.path(out_dir, "all_sample_h1_h3_pathways.csv"), row.names = FALSE)
write.csv(h_summary, file.path(out_dir, "h1_vs_h3_pathway_summary.csv"), row.names = FALSE)
write.csv(c, file.path(out_dir, "all_sample_competitive_h3.csv"), row.names = FALSE)
write.csv(c_summary, file.path(out_dir, "competitive_h3_pathway_summary.csv"), row.names = FALSE)

md <- c(
  "# PASSAGE Cancer Panel Summary",
  "",
  paste0("- PASSAGE input rows: ", nrow(manifest)),
  paste0("- H1/H3 files found: ", length(h_files[file.exists(h_files)]), " / ", nrow(manifest)),
  paste0("- Competitive H3 files found: ", length(c_files[file.exists(c_files)]), " / ", nrow(manifest)),
  "",
  "Primary H3 covariates are technical covariates plus broad marker-derived cell-type scores.",
  "Raw scRNA references are staged under refs/raw_scRNA for the next deconvolution pass.",
  "Competitive H3 uses matched-module Monte Carlo p-values for score_z and score_robust_z; cEPSV and PC1 coherence are reported as descriptors only.",
  "",
  "## Top Competitive H3",
  ""
)
if (nrow(c_summary)) {
  top <- head(c_summary[, c("cancer", "reference", "statistic", "pathway",
                            "median_p", "fisher_p", "fisher_fdr",
                            "median_cEPSV", "median_coherence_pc1")], 30)
  md <- c(md, paste(colnames(top), collapse = " | "),
          paste(rep("---", ncol(top)), collapse = " | "),
          apply(top, 1L, function(z) paste(z, collapse = " | ")))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote aggregate outputs to ", out_dir)
