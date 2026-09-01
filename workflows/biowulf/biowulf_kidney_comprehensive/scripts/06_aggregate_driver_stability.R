#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
res_root <- file.path(root, "results", "driver_stability")
out_dir <- file.path(root, "results", "driver_stability_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_many <- function(pattern) {
  files <- list.files(res_root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(data.frame())
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

summ <- read_many("_summary[.]csv$")
drivers <- read_many("_drivers[.]csv$")
nulls <- read_many("_nulls[.]csv$")
write.csv(summ, file.path(out_dir, "driver_stability_all_summaries.csv"), row.names = FALSE)
write.csv(drivers, file.path(out_dir, "driver_stability_all_driver_bootstraps.csv"), row.names = FALSE)
write.csv(nulls, file.path(out_dir, "driver_stability_all_nulls.csv"), row.names = FALSE)

if (nrow(summ)) {
  ref_stat <- do.call(rbind, lapply(split(summ, list(summ$reference, summ$statistic), drop = TRUE), function(z) {
    data.frame(
      reference = z$reference[[1L]],
      statistic = z$statistic[[1L]],
      n_pathway_samples = nrow(z),
      median_max_selection_frequency = stats::median(z$max_selection_frequency, na.rm = TRUE),
      median_topk_jaccard = stats::median(z$mean_topk_jaccard, na.rm = TRUE),
      median_n_genes_freq_ge_050 = stats::median(z$n_genes_freq_ge_050, na.rm = TRUE),
      frac_null_p_maxfreq_le_020 = mean(z$null_p_max_frequency <= 0.20, na.rm = TRUE),
      frac_null_p_jaccard_le_020 = mean(z$null_p_jaccard <= 0.20, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(ref_stat, file.path(out_dir, "driver_stability_by_reference_statistic.csv"), row.names = FALSE)
} else {
  ref_stat <- data.frame()
}

if (nrow(drivers)) {
  gene_key <- interaction(drivers$sample, drivers$reference, drivers$pathway, drivers$statistic, drivers$gene, drop = TRUE, sep = "\r")
  gene_stab <- do.call(rbind, lapply(split(drivers, gene_key), function(z) {
    data.frame(
      sample = z$sample[[1L]],
      reference = z$reference[[1L]],
      pathway = z$pathway[[1L]],
      statistic = z$statistic[[1L]],
      gene = z$gene[[1L]],
      selection_frequency = length(unique(z$bootstrap)) / max(z$bootstrap, na.rm = TRUE),
      mean_rank = mean(z$rank, na.rm = TRUE),
      median_weight = stats::median(z$driver_weight, na.rm = TRUE),
      mean_raw_score = mean(z$raw_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  gene_stab <- gene_stab[order(gene_stab$reference, gene_stab$sample, gene_stab$pathway,
                               gene_stab$statistic, -gene_stab$selection_frequency, gene_stab$mean_rank), , drop = FALSE]
  write.csv(gene_stab, file.path(out_dir, "driver_stability_gene_frequencies.csv"), row.names = FALSE)
} else {
  gene_stab <- data.frame()
}

top_stable <- if (nrow(gene_stab)) {
  z <- gene_stab[gene_stab$selection_frequency >= 0.50, , drop = FALSE]
  z <- z[order(z$reference, z$statistic, z$pathway, -z$selection_frequency, z$mean_rank), , drop = FALSE]
  z
} else data.frame()
write.csv(top_stable, file.path(out_dir, "driver_stability_stable_genes_freq_ge_050.csv"), row.names = FALSE)

md <- c(
  "# Kidney PASSAGE Driver-Gene Stability",
  "",
  paste0("- summary rows: ", nrow(summ)),
  paste0("- driver bootstrap rows: ", nrow(drivers)),
  paste0("- null rows: ", nrow(nulls)),
  paste0("- stable gene rows freq>=0.50: ", nrow(top_stable)),
  "",
  "## Reference/Statistic Summary",
  "",
  paste(c("reference", "statistic", "n", "median_max_freq", "median_jaccard", "median_n_freq_ge_050", "frac_null_p_maxfreq_le_020", "frac_null_p_jaccard_le_020"), collapse = " | "),
  paste(rep("---", 8), collapse = " | ")
)
if (nrow(ref_stat)) {
  for (ii in seq_len(nrow(ref_stat))) {
    r <- ref_stat[ii, ]
    md <- c(md, paste(c(r$reference, r$statistic, r$n_pathway_samples,
                        sprintf("%.3f", r$median_max_selection_frequency),
                        sprintf("%.3f", r$median_topk_jaccard),
                        sprintf("%.1f", r$median_n_genes_freq_ge_050),
                        sprintf("%.3f", r$frac_null_p_maxfreq_le_020),
                        sprintf("%.3f", r$frac_null_p_jaccard_le_020)), collapse = " | "))
  }
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
