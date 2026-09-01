#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)

summary_dir <- file.path(root, "results", "passage_kidney_summary")
out_dir <- file.path(root, "results", "passage_gwas_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comp_file <- file.path(summary_dir, "all_sample_competitive_h3.csv")
h_file <- file.path(summary_dir, "all_sample_h1_h3_pathways.csv")
if (!file.exists(comp_file)) stop("Missing ", comp_file)
if (!file.exists(h_file)) stop("Missing ", h_file)

comp <- read.csv(comp_file, stringsAsFactors = FALSE)
h <- read.csv(h_file, stringsAsFactors = FALSE)
scorez <- comp[comp$statistic == "score_z", , drop = FALSE]

split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  unique(unlist(strsplit(x[nzchar(x)], ";", fixed = TRUE), use.names = FALSE))
}

driver_rows <- list()
for (ii in seq_len(nrow(h))) {
  genes <- split_genes(h$spasset_genes_H3[[ii]])
  if (!length(genes)) next
  driver_rows[[length(driver_rows) + 1L]] <- data.frame(
    sample = h$sample[[ii]],
    spatial_sample = h$spatial_sample[[ii]],
    reference = h$reference[[ii]],
    pathway = h$pathway[[ii]],
    gene = genes,
    p_spasset_H3 = h$p_spasset_H3[[ii]],
    spasset_size_H3 = h$spasset_size_H3[[ii]],
    p_H3 = h$p_H3[[ii]],
    fdr_H3 = h$fdr_H3[[ii]],
    stringsAsFactors = FALSE
  )
}
drivers <- if (length(driver_rows)) do.call(rbind, driver_rows) else data.frame()
write.csv(drivers, file.path(out_dir, "passage_spasset_driver_gene_long.csv"), row.names = FALSE)

priority <- scorez
priority$primary_passage <- priority$competitive_score_fdr_global <= 0.05
priority$empirical_passage <- priority$competitive_score_empirical_leave_sample_fdr <= 0.05
priority <- priority[order(priority$reference, priority$competitive_score_p), , drop = FALSE]
write.csv(priority, file.path(out_dir, "passage_scorez_pathway_priority_table.csv"), row.names = FALSE)

write_gmt <- function(named_sets, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(named_sets)) {
    genes <- sort(unique(named_sets[[nm]]))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    if (length(genes) >= 2L) {
      writeLines(paste(c(nm, "PASSAGE", genes), collapse = "\t"), con)
    }
  }
}

driver_sets <- list()
if (nrow(drivers)) {
  key <- paste(drivers$reference, drivers$pathway, sep = "__")
  driver_sets <- split(drivers$gene, key)
}
write_gmt(driver_sets, file.path(out_dir, "passage_spasset_driver_sets.gmt"))

sig <- priority[priority$competitive_score_fdr_global <= 0.05, , drop = FALSE]
sig_driver_sets <- list()
if (nrow(sig) && nrow(drivers)) {
  sig_keys <- paste(sig$reference, sig$pathway, sep = "__")
  keep <- paste(drivers$reference, drivers$pathway, sep = "__") %in% sig_keys
  sig_driver_sets <- split(drivers$gene[keep], paste(drivers$reference[keep], drivers$pathway[keep], sep = "__"))
}
write_gmt(sig_driver_sets, file.path(out_dir, "passage_significant_scorez_driver_sets.gmt"))

recurrent <- do.call(rbind, lapply(split(priority, paste(priority$reference, priority$pathway, sep = "__")), function(z) {
  data.frame(
    reference = z$reference[[1L]],
    pathway = z$pathway[[1L]],
    n_samples = length(unique(z$spatial_sample)),
    n_raw_p05 = sum(z$competitive_score_p <= 0.05, na.rm = TRUE),
    n_fdr05 = sum(z$competitive_score_fdr_global <= 0.05, na.rm = TRUE),
    n_empirical_fdr05 = sum(z$competitive_score_empirical_leave_sample_fdr <= 0.05, na.rm = TRUE),
    median_p = stats::median(z$competitive_score_p, na.rm = TRUE),
    median_fdr = stats::median(z$competitive_score_fdr_global, na.rm = TRUE),
    median_cEPSV = stats::median(z$cEPSV, na.rm = TRUE),
    median_coherence_pc1 = stats::median(z$coherence_pc1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
recurrent <- recurrent[order(recurrent$reference, recurrent$median_p), , drop = FALSE]
write.csv(recurrent, file.path(out_dir, "passage_recurrent_pathway_summary_for_gwas.csv"), row.names = FALSE)

manifest <- file.path(root, "refs", "gwas", "rcc_gwas_sumstats_manifest.csv")
if (file.exists(manifest)) {
  file.copy(manifest, file.path(out_dir, "rcc_gwas_sumstats_manifest.csv"), overwrite = TRUE)
}

writeLines(c(
  "# PASSAGE-GWAS Inputs",
  "",
  paste0("- driver_long_rows: ", nrow(drivers)),
  paste0("- scorez_rows: ", nrow(priority)),
  paste0("- significant_scorez_rows: ", sum(priority$competitive_score_fdr_global <= 0.05, na.rm = TRUE)),
  "",
  "Files:",
  "- passage_spasset_driver_gene_long.csv",
  "- passage_scorez_pathway_priority_table.csv",
  "- passage_spasset_driver_sets.gmt",
  "- passage_significant_scorez_driver_sets.gmt",
  "- passage_recurrent_pathway_summary_for_gwas.csv",
  "- rcc_gwas_sumstats_manifest.csv"
), file.path(out_dir, "README.md"))

message("Wrote PASSAGE-GWAS input files to ", out_dir)
