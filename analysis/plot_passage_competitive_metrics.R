# Summarize and render plots for PASSAGE competitive metric outputs.

parse_args <- function(args) {
  cfg <- list(
    result_dirs = character(),
    labels = character(),
    condition = "H3",
    out_dir = "results/passage_competitive_plots"
  )
  for (arg in args) {
    if (grepl("^--result-dirs=", arg)) cfg$result_dirs <- strsplit(sub("^--result-dirs=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--labels=", arg)) cfg$labels <- strsplit(sub("^--labels=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--condition=", arg)) cfg$condition <- sub("^--condition=", "", arg)
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
  }
  if (length(cfg$result_dirs) == 0L) stop("--result-dirs is required")
  if (length(cfg$labels) == 0L) {
    cfg$labels <- basename(dirname(cfg$result_dirs))
  }
  if (length(cfg$labels) != length(cfg$result_dirs)) stop("--labels must align with --result-dirs")
  cfg
}

safe_p <- function(x) pmax(as.numeric(x), .Machine$double.xmin)
neglog10 <- function(x) -log10(safe_p(x))

read_one <- function(dir, label, condition) {
  comp_path <- file.path(dir, paste0("passage_competitive_metrics_", condition, ".csv"))
  if (!file.exists(comp_path)) stop("Missing competitive metrics: ", comp_path)
  comp <- read.csv(comp_path, stringsAsFactors = FALSE)
  comp$analysis_label <- label
  base_path <- file.path(dir, "passage_msigdb_conditional_pathways.csv")
  if (file.exists(base_path)) {
    base <- read.csv(base_path, stringsAsFactors = FALSE)
    keep <- intersect(c("pathway", "p_H1", "p_H2", "p_H3", "p_H4",
                        "fdr_H1", "fdr_H2", "fdr_H3", "fdr_H4",
                        "svg_acat_p", "svg_acat_fdr",
                        "cell_type_share", "background_share", "pathway_specific_share"),
                      names(base))
    comp <- merge(comp, base[, keep, drop = FALSE], by = "pathway", all.x = TRUE)
  }
  comp
}

count_sig <- function(x, col) {
  if (!col %in% names(x)) return(NA_integer_)
  sum(is.finite(x[[col]]) & x[[col]] < 0.05)
}

plot_counts <- function(tbl, out_dir, condition) {
  labs <- unique(tbl$analysis_label)
  base_col <- paste0("fdr_", condition)
  mat <- sapply(labs, function(lb) {
    x <- tbl[tbl$analysis_label == lb, , drop = FALSE]
    c(
      PASSAGE = count_sig(x, base_col),
      meanPVE_z = count_sig(x, "competitive_mean_propSV_fdr_analytic"),
      meanPVE_perm = count_sig(x, "competitive_mean_propSV_fdr_perm"),
      cEPSV_perm = count_sig(x, "competitive_cEPSV_fdr_perm"),
      coherence_perm = count_sig(x, "coherence_pc1_fdr_perm")
    )
  })
  png(file.path(out_dir, "significant_counts.png"), width = 1500, height = 850, res = 140)
  par(mar = c(9, 5, 3, 1))
  barplot(mat, beside = TRUE, las = 2, col = c("#4C78A8", "#F58518", "#E45756", "#54A24B", "#B279A2"),
          ylab = "FDR < 0.05 pathways", main = paste(condition, "PASSAGE and Competitive Null Discoveries"))
  legend("topright", legend = rownames(mat), fill = c("#4C78A8", "#F58518", "#E45756", "#54A24B", "#B279A2"), cex = 0.8)
  dev.off()
}

plot_scatter <- function(tbl, out_dir, condition) {
  p_col <- paste0("p_", condition)
  for (lb in unique(tbl$analysis_label)) {
    x <- tbl[tbl$analysis_label == lb, , drop = FALSE]
    slug <- gsub("[^A-Za-z0-9]+", "_", lb)
    if (all(c(p_col, "competitive_cEPSV_p_perm") %in% names(x))) {
      png(file.path(out_dir, paste0("scatter_", tolower(condition), "_vs_cepsv_", slug, ".png")), width = 950, height = 850, res = 130)
      sig <- is.finite(x$competitive_cEPSV_fdr_perm) & x$competitive_cEPSV_fdr_perm < 0.05
      plot(neglog10(x[[p_col]]), neglog10(x$competitive_cEPSV_p_perm),
           pch = 19, col = ifelse(sig, "#D62728AA", "#2F4B7CAA"),
           xlab = paste0("-log10 PASSAGE ", condition, " p-value"), ylab = "-log10 competitive cEPSV permutation p-value",
           main = lb)
      abline(h = -log10(0.05), v = -log10(0.05), lty = 2, col = "gray50")
      legend("topleft", legend = c("competitive cEPSV FDR < 0.05", "not significant"),
             col = c("#D62728AA", "#2F4B7CAA"), pch = 19, bty = "n", cex = 0.85)
      dev.off()
    }
    if (all(c("pathway_specific_share", "cEPSV") %in% names(x))) {
      png(file.path(out_dir, paste0("scatter_specific_share_cepsv_", slug, ".png")), width = 950, height = 850, res = 130)
      sig <- is.finite(x$competitive_cEPSV_fdr_perm) & x$competitive_cEPSV_fdr_perm < 0.05
      plot(x$pathway_specific_share, x$cEPSV,
           pch = 19, col = ifelse(sig, "#D62728AA", "#2F4B7CAA"),
           xlab = "PASSAGE pathway-specific share (Q_H3 / Q_H1)",
           ylab = "cEPSV", main = lb)
      legend("topright", legend = c("competitive cEPSV FDR < 0.05", "not significant"),
             col = c("#D62728AA", "#2F4B7CAA"), pch = 19, bty = "n", cex = 0.85)
      dev.off()
    }
  }
}

plot_top <- function(tbl, out_dir, n = 20L) {
  for (lb in unique(tbl$analysis_label)) {
    x <- tbl[tbl$analysis_label == lb & is.finite(tbl$cEPSV), , drop = FALSE]
    if (nrow(x) == 0L) next
    x <- x[order(x$competitive_cEPSV_p_perm, -x$cEPSV, na.last = TRUE), , drop = FALSE]
    x <- head(x, n)
    slug <- gsub("[^A-Za-z0-9]+", "_", lb)
    labs <- gsub("^HALLMARK_|^KEGG_", "", x$pathway)
    labs <- substr(labs, 1, 55)
    png(file.path(out_dir, paste0("top_cepsv_", slug, ".png")), width = 1300, height = 950, res = 130)
    par(mar = c(5, 18, 3, 1))
    barplot(rev(x$cEPSV), horiz = TRUE, names.arg = rev(labs), las = 1,
            col = ifelse(rev(x$competitive_cEPSV_fdr_perm) < 0.05, "#D62728", "#4C78A8"),
            xlab = "cEPSV", main = paste("Top competitive cEPSV pathways:", lb))
    dev.off()
  }
}

write_summary <- function(tbl, out_dir, condition) {
  fdr_col <- paste0("fdr_", condition)
  lines <- c("# PASSAGE Competitive Metrics Summary", "")
  for (lb in unique(tbl$analysis_label)) {
    x <- tbl[tbl$analysis_label == lb, , drop = FALSE]
    lines <- c(lines,
      paste0("## ", lb),
      paste0("- Pathways: ", nrow(x)),
      paste0("- PASSAGE ", condition, " FDR < 0.05: ", count_sig(x, fdr_col)),
      paste0("- Competitive mean propSV analytic FDR < 0.05: ", count_sig(x, "competitive_mean_propSV_fdr_analytic")),
      paste0("- Competitive mean propSV permutation FDR < 0.05: ", count_sig(x, "competitive_mean_propSV_fdr_perm")),
      paste0("- Competitive cEPSV permutation FDR < 0.05: ", count_sig(x, "competitive_cEPSV_fdr_perm")),
      paste0("- Pathway coherence permutation FDR < 0.05: ", count_sig(x, "coherence_pc1_fdr_perm")),
      ""
    )
    top <- head(x[order(x$competitive_cEPSV_p_perm, -x$cEPSV, na.last = TRUE),
                  c("pathway", "cEPSV", "competitive_cEPSV_p_perm", "competitive_cEPSV_fdr_perm",
                    "pc1_spatial_fraction", "coherence_pc1_p_perm", "coherence_pc1_fdr_perm"),
                  drop = FALSE], 10)
    lines <- c(lines, "pathway | cEPSV | cEPSV p | cEPSV FDR | PC1 spatial fraction | coherence p | coherence FDR",
               "--- | ---: | ---: | ---: | ---: | ---: | ---:")
    if (nrow(top) > 0L) {
      for (ii in seq_len(nrow(top))) {
        lines <- c(lines, paste(top[ii, ], collapse = " | "))
      }
    }
    lines <- c(lines, "")
  }
  writeLines(lines, file.path(out_dir, "competitive_metrics_summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
tbl <- do.call(rbind, Map(read_one, cfg$result_dirs, cfg$labels, MoreArgs = list(condition = cfg$condition)))
write.csv(tbl, file.path(cfg$out_dir, "competitive_metrics_combined.csv"), row.names = FALSE)
plot_counts(tbl, cfg$out_dir, cfg$condition)
plot_scatter(tbl, cfg$out_dir, cfg$condition)
plot_top(tbl, cfg$out_dir)
write_summary(tbl, cfg$out_dir, cfg$condition)
message("Wrote plots and summary to ", cfg$out_dir)
