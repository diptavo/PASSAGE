#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
top_k <- if (length(args) >= 2L) as.integer(args[[2L]]) else 10L

sum_dir <- file.path(root, "results", "driver_stability_summary")
gene_file <- file.path(sum_dir, "driver_stability_gene_frequencies.csv")
summ_file <- file.path(sum_dir, "driver_stability_all_summaries.csv")
if (!file.exists(gene_file)) stop("Missing gene frequency table: ", gene_file)
if (!file.exists(summ_file)) stop("Missing summary table: ", summ_file)

gene <- read.csv(gene_file, stringsAsFactors = FALSE)
summ <- read.csv(summ_file, stringsAsFactors = FALSE)
stat_order <- c("CSPS", "GSPS", "MMP", "CSV", "HCPS", "OTSAS", "score_z", "score_z_robust")

split_id <- interaction(gene$cohort, gene$sample, gene$pathway, gene$statistic, drop = TRUE, sep = "\r")
top_sets <- lapply(split(gene, split_id), function(z) {
  z <- z[order(-z$selection_frequency, z$mean_rank), , drop = FALSE]
  head(z$gene, top_k)
})
parse_key <- function(x) strsplit(as.character(x), "\r", fixed = TRUE)[[1L]]

pairs <- list()
ii <- 0L
groups <- unique(interaction(summ$cohort, summ$sample, summ$pathway, drop = TRUE, sep = "\r"))
for (g in groups) {
  parts <- parse_key(g)
  cohort <- parts[[1L]]
  sample <- parts[[2L]]
  pathway <- parts[[3L]]
  keys <- paste(cohort, sample, pathway, stat_order, sep = "\r")
  present <- stat_order[keys %in% names(top_sets)]
  if (length(present) < 2L) next
  for (a in seq_len(length(present) - 1L)) {
    for (b in seq.int(a + 1L, length(present))) {
      s1 <- present[[a]]
      s2 <- present[[b]]
      set1 <- top_sets[[paste(cohort, sample, pathway, s1, sep = "\r")]]
      set2 <- top_sets[[paste(cohort, sample, pathway, s2, sep = "\r")]]
      u <- union(set1, set2)
      ii <- ii + 1L
      pairs[[ii]] <- data.frame(
        cohort = cohort,
        sample = sample,
        pathway = pathway,
        stat1 = s1,
        stat2 = s2,
        top_k = top_k,
        overlap = length(intersect(set1, set2)),
        jaccard = if (length(u)) length(intersect(set1, set2)) / length(u) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
pair_tbl <- do.call(rbind, pairs)
write.csv(pair_tbl, file.path(sum_dir, "driver_stat_pairwise_topk_overlap.csv"), row.names = FALSE)

pair_sum <- do.call(rbind, lapply(split(pair_tbl, list(pair_tbl$cohort, pair_tbl$stat1, pair_tbl$stat2), drop = TRUE), function(z) {
  data.frame(
    cohort = z$cohort[[1L]],
    stat1 = z$stat1[[1L]],
    stat2 = z$stat2[[1L]],
    n = nrow(z),
    median_overlap = stats::median(z$overlap, na.rm = TRUE),
    median_jaccard = stats::median(z$jaccard, na.rm = TRUE),
    mean_jaccard = mean(z$jaccard, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
pair_sum <- pair_sum[order(pair_sum$cohort, -pair_sum$median_jaccard, pair_sum$stat1, pair_sum$stat2), , drop = FALSE]
write.csv(pair_sum, file.path(sum_dir, "driver_stat_pairwise_topk_overlap_summary.csv"), row.names = FALSE)

baseline <- pair_sum[pair_sum$stat1 %in% c("CSPS", "GSPS", "MMP", "CSV", "HCPS", "OTSAS") &
                       pair_sum$stat2 %in% c("score_z", "score_z_robust"), , drop = FALSE]
write.csv(baseline, file.path(sum_dir, "driver_stat_overlap_vs_scorez_summary.csv"), row.names = FALSE)

md <- c(
  "# Driver Agreement Across Statistics",
  "",
  paste0("- top_k: ", top_k),
  paste0("- pairwise rows: ", nrow(pair_tbl)),
  "",
  "## Highest Median Jaccard Per Cohort",
  "",
  paste(c("cohort", "stat1", "stat2", "n", "median_overlap", "median_jaccard"), collapse = " | "),
  paste(rep("---", 6), collapse = " | ")
)
for (cc in sort(unique(pair_sum$cohort))) {
  z <- head(pair_sum[pair_sum$cohort == cc, , drop = FALSE], 12L)
  for (rr in seq_len(nrow(z))) {
    r <- z[rr, ]
    md <- c(md, paste(c(r$cohort, r$stat1, r$stat2, r$n,
                        sprintf("%.1f", r$median_overlap),
                        sprintf("%.3f", r$median_jaccard)), collapse = " | "))
  }
}
writeLines(md, file.path(sum_dir, "driver_stat_agreement.md"))
message("Wrote driver statistic agreement outputs to ", sum_dir)
