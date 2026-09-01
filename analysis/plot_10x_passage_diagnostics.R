# Plot diagnostics for PASSAGE Hallmark runs on the 10x Visium breast cancer data.
#
# Usage:
#   Rscript scripts/plot_10x_passage_diagnostics.R \
#     --result-root=results/passage_10x_hallmark_perm999_fastcal

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
})

parse_args <- function(args) {
  cfg <- list(
    result_root = "results/passage_10x_hallmark_perm999_fastcal",
    top_n_pathways = 8L,
    max_driver_genes = 12L
  )
  for (arg in args) {
    if (grepl("^--result-root=", arg)) cfg$result_root <- sub("^--result-root=", "", arg)
    if (grepl("^--top-n-pathways=", arg)) cfg$top_n_pathways <- as.integer(sub("^--top-n-pathways=", "", arg))
    if (grepl("^--max-driver-genes=", arg)) cfg$max_driver_genes <- as.integer(sub("^--max-driver-genes=", "", arg))
  }
  cfg
}

plot_spatial_value <- function(df, value_col, title, out_file, point_size = 0.55) {
  p <- ggplot(df, aes(x = x, y = y, color = .data[[value_col]])) +
    geom_point(size = point_size, alpha = 0.9) +
    scale_color_viridis(option = "magma") +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, color = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      legend.position = "right"
    )
  ggsave(out_file, p, width = 5.2, height = 4.6, dpi = 180)
}

zscore_cols <- function(x) {
  x <- as.matrix(x)
  s <- apply(x, 2L, stats::sd)
  s <- pmax(s, sqrt(.Machine$double.eps))
  sweep(sweep(x, 2L, colMeans(x), "-"), 2L, s, "/")
}

pathway_score <- function(Y, genes) {
  genes <- intersect(genes, colnames(Y))
  if (length(genes) < 1L) {
    return(rep(NA_real_, nrow(Y)))
  }
  rowMeans(zscore_cols(Y[, genes, drop = FALSE]))
}

make_base_df <- function(dat) {
  md <- dat$spot_metadata
  x_col <- if ("pxl_col_in_fullres" %in% names(md)) "pxl_col_in_fullres" else "array_col"
  y_col <- if ("pxl_row_in_fullres" %in% names(md)) "pxl_row_in_fullres" else "array_row"
  data.frame(
    barcode = md$barcode,
    x = md[[x_col]],
    y = md[[y_col]],
    log_umi = log1p(md$lib_size),
    log_detected = log1p(md$n_detected),
    array_row = md$array_row,
    array_col = md$array_col,
    stringsAsFactors = FALSE
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

summarize_correlations <- function(df, variables, covariates) {
  rows <- list()
  for (v in variables) {
    for (cvar in covariates) {
      ok <- is.finite(df[[v]]) & is.finite(df[[cvar]])
      cc <- if (sum(ok) >= 3L) stats::cor(df[[v]][ok], df[[cvar]][ok]) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        variable = v,
        covariate = cvar,
        correlation = cc,
        abs_correlation = abs(cc),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_one_dataset <- function(dataset_dir, cfg) {
  dataset_id <- basename(dataset_dir)
  message("Diagnostics for ", dataset_id)

  fit <- readRDS(file.path(dataset_dir, "passage_hallmark_result.rds"))
  dat <- readRDS(file.path(dataset_dir, "passage_hallmark_prepared_data.rds"))
  out_dir <- file.path(dataset_dir, "diagnostics")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  base <- make_base_df(dat)
  plot_spatial_value(base, "log_umi", paste(dataset_id, "log UMI"), file.path(out_dir, "qc_log_umi.png"))
  plot_spatial_value(base, "log_detected", paste(dataset_id, "log detected genes"), file.path(out_dir, "qc_log_detected.png"))

  v <- as.data.frame(fit$engine$V)
  names(v) <- paste0("factor_", seq_len(ncol(v)))
  df <- cbind(base, v)
  for (k in seq_len(min(6L, ncol(v)))) {
    plot_spatial_value(
      df,
      paste0("factor_", k),
      paste(dataset_id, "PASSAGE factor", k),
      file.path(out_dir, paste0("factor_", k, ".png"))
    )
  }

  tbl <- fit$summary
  tbl <- tbl[order(tbl$p_H1, -tbl$R2_cca, tbl$pathway), , drop = FALSE]
  top_pathways <- head(tbl$pathway, cfg$top_n_pathways)
  score_names <- character()
  for (pname in top_pathways) {
    score_col <- paste0("score_", safe_name(pname))
    df[[score_col]] <- pathway_score(dat$Y, dat$pathways[[pname]])
    score_names <- c(score_names, score_col)
    plot_spatial_value(
      df,
      score_col,
      paste(dataset_id, pname),
      file.path(out_dir, paste0("pathway_", safe_name(pname), ".png"))
    )
  }

  driver_genes <- unique(unlist(strsplit(head(tbl$spasset_genes, cfg$top_n_pathways), ";", fixed = TRUE)))
  driver_genes <- driver_genes[nzchar(driver_genes)]
  driver_genes <- head(intersect(driver_genes, colnames(dat$Y)), cfg$max_driver_genes)
  for (gene in driver_genes) {
    gene_col <- paste0("gene_", safe_name(gene))
    df[[gene_col]] <- dat$Y[, gene]
    plot_spatial_value(
      df,
      gene_col,
      paste(dataset_id, gene),
      file.path(out_dir, paste0("driver_gene_", safe_name(gene), ".png"))
    )
  }

  covars <- c("log_umi", "log_detected", "x", "y", "array_row", "array_col")
  variables <- c(names(v), score_names)
  cor_tbl <- summarize_correlations(df, variables, covars)
  cor_tbl <- cor_tbl[order(-cor_tbl$abs_correlation), , drop = FALSE]
  write.csv(cor_tbl, file.path(out_dir, "diagnostic_correlations.csv"), row.names = FALSE)

  top_cor <- head(cor_tbl, 12)
  lines <- c(
    paste0("# Diagnostics - ", dataset_id),
    "",
    paste0("- Spots: ", nrow(dat$Y)),
    paste0("- Hallmark genes: ", ncol(dat$Y)),
    paste0("- Hallmark pathways: ", nrow(tbl)),
    paste0("- FDR <= 0.05 pathways: ", sum(tbl$fdr_H1 <= 0.05)),
    "",
    "## Plotted Pathways",
    "",
    paste0("- ", top_pathways),
    "",
    "## Driver Genes Plotted",
    "",
    if (length(driver_genes) > 0L) paste0("- ", driver_genes) else "- none",
    "",
    "## Largest Absolute Diagnostic Correlations",
    "",
    "variable | covariate | correlation",
    "--- | --- | ---"
  )
  cor_lines <- apply(top_cor[, c("variable", "covariate", "correlation")], 1L, function(r) {
    paste(r, collapse = " | ")
  })
  writeLines(c(lines, cor_lines), file.path(out_dir, "diagnostic_summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dataset_dirs <- list.dirs(cfg$result_root, recursive = FALSE, full.names = TRUE)
if (length(dataset_dirs) == 0L) {
  stop("No dataset result folders found under ", cfg$result_root)
}
for (dataset_dir in dataset_dirs) {
  run_one_dataset(dataset_dir, cfg)
}
