#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)
for (f in sort(list.files("SpaPath/R", pattern = "[.]R$", full.names = TRUE))) source(f)

out_dir <- file.path(root, "results", "passage_dlpfc_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_many <- function(files) {
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

manifest <- read.csv(file.path(root, "data", "passage_inputs", "passage_input_manifest.csv"),
                     stringsAsFactors = FALSE)
h_files <- file.path(root, "results", "passage_dlpfc_h1_h3", manifest$sample,
                     "passage_dlpfc_h1_h3_pathways.csv")
c_files <- file.path(root, "results", "passage_dlpfc_competitive_h3", manifest$sample,
                     "passage_dlpfc_competitive_h3.csv")
h <- read_many(h_files)
c <- read_many(c_files)

combine_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}
summarize_h <- function(x) {
  data.frame(
    n_sections = length(unique(x$sample)),
    n_H1_fdr05 = sum(x$fdr_H1 <= 0.05, na.rm = TRUE),
    n_H3_fdr05 = sum(x$fdr_H3 <= 0.05, na.rm = TRUE),
    median_p_H1 = stats::median(x$p_H1, na.rm = TRUE),
    median_p_H3 = stats::median(x$p_H3, na.rm = TRUE),
    fisher_p_H1 = combine_p(x$p_H1),
    fisher_p_H3 = combine_p(x$p_H3),
    median_celltype_adjusted_reduction = stats::median(x$celltype_adjusted_reduction, na.rm = TRUE),
    median_H3_over_H1_Q = stats::median(x$H3_over_H1_Q, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
h_summary <- if (nrow(h)) {
  do.call(rbind, lapply(split(h, h$pathway), function(z) {
    cbind(pathway = z$pathway[[1L]], summarize_h(z), stringsAsFactors = FALSE)
  }))
} else data.frame()
if (nrow(h_summary)) {
  h_summary$fisher_fdr_H1 <- stats::p.adjust(h_summary$fisher_p_H1, method = "BH")
  h_summary$fisher_fdr_H3 <- stats::p.adjust(h_summary$fisher_p_H3, method = "BH")
  h_summary <- h_summary[order(h_summary$fisher_p_H3, h_summary$fisher_p_H1), , drop = FALSE]
}

if (nrow(c)) {
  c$dataset <- "spatialDLPFC"
  c$mode <- "technical_engine_h3_score"
  c$matching <- "expr_detect_var_spatial"
  c$mc_sampler <- "module"
  c$replicate <- c$sample
  c <- passage_empirical_tail_calibrate(
    c, p_col = "competitive_score_p",
    group_cols = c("dataset", "statistic", "mode", "matching", "mc_sampler"),
    leaveout_cols = "subject",
    out_col = "competitive_score_empirical_leave_subject_p",
    n_col = "competitive_score_empirical_leave_subject_n"
  )
  c$competitive_score_fdr_global <- ave(c$competitive_score_p, c$statistic,
                                        FUN = function(p) stats::p.adjust(p, method = "BH"))
  c$competitive_score_empirical_leave_subject_fdr <- ave(
    c$competitive_score_empirical_leave_subject_p, c$statistic,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
}

write.csv(h, file.path(out_dir, "all_section_h1_h3_pathways.csv"), row.names = FALSE)
write.csv(h_summary, file.path(out_dir, "h1_vs_h3_pathway_summary.csv"), row.names = FALSE)
write.csv(c, file.path(out_dir, "all_section_competitive_h3.csv"), row.names = FALSE)

md <- c(
  "# PASSAGE spatialDLPFC Summary",
  "",
  paste0("- H1/H3 section files found: ", length(h_files[file.exists(h_files)]), " / ", nrow(manifest)),
  paste0("- Competitive section files found: ", length(c_files[file.exists(c_files)]), " / ", nrow(manifest)),
  "",
  "Primary H3 covariates are technical covariates plus broad cell2location proportions.",
  "Competitive H3 uses raw matched-module p-values plus a leave-subject empirical ranking column.",
  "",
  "## Top H3 Fisher Pathways",
  ""
)
if (nrow(h_summary)) {
  top <- head(h_summary[, c("pathway", "n_H1_fdr05", "n_H3_fdr05", "fisher_p_H1",
                            "fisher_p_H3", "fisher_fdr_H3",
                            "median_celltype_adjusted_reduction")], 20)
  md <- c(md, paste(colnames(top), collapse = " | "),
          paste(rep("---", ncol(top)), collapse = " | "),
          apply(top, 1L, function(z) paste(z, collapse = " | ")))
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote aggregate outputs to ", out_dir)
