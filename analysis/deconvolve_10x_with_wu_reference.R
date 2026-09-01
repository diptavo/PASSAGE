# Estimate spot-level cell-type proportions for 10x Visium breast cancer data
# using the Wu et al. GSE176078 coarse pseudobulk reference.
#
# This uses a fast projected-gradient fit:
#   min_W || W S - Y ||^2, subject to W rows on the simplex.
#
# It also saves broad background factors from non-Hallmark variable genes for
# downstream H3 conditioning.

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(msigdbr)
  library(ggplot2)
  library(viridis)
})

parse_args <- function(args) {
  cfg <- list(
    reference_rds = "data/reference/GSE176078/processed/wu_brca_coarse_reference.rds",
    out_root = "results/wu_deconvolution_10x",
    datasets = c("Visium_FFPE_Human_Breast_Cancer", "V1_Breast_Cancer_Block_A_Section_1"),
    marker_per_type = 150L,
    min_marker_score = 0.10,
    n_bg_factors = 6L,
    n_bg_genes = 2500L,
    deconv_iter = 350L
  )
  for (arg in args) {
    if (grepl("^--reference-rds=", arg)) cfg$reference_rds <- sub("^--reference-rds=", "", arg)
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--marker-per-type=", arg)) cfg$marker_per_type <- as.integer(sub("^--marker-per-type=", "", arg))
    if (grepl("^--min-marker-score=", arg)) cfg$min_marker_score <- as.numeric(sub("^--min-marker-score=", "", arg))
    if (grepl("^--n-bg-factors=", arg)) cfg$n_bg_factors <- as.integer(sub("^--n-bg-factors=", "", arg))
    if (grepl("^--n-bg-genes=", arg)) cfg$n_bg_genes <- as.integer(sub("^--n-bg-genes=", "", arg))
    if (grepl("^--deconv-iter=", arg)) cfg$deconv_iter <- as.integer(sub("^--deconv-iter=", "", arg))
  }
  cfg
}

dataset_configs <- list(
  Visium_FFPE_Human_Breast_Cancer = list(
    label = "10x Visium FFPE Human Breast Cancer DCIS/Invasive",
    root = file.path("data", "raw", "Visium_FFPE_Human_Breast_Cancer")
  ),
  V1_Breast_Cancer_Block_A_Section_1 = list(
    label = "10x Visium Breast Cancer Block A Section 1",
    root = file.path("data", "raw", "V1_Breast_Cancer_Block_A_Section_1")
  )
)

read_10x_h5 <- function(h5_path) {
  data <- rhdf5::h5read(h5_path, "/matrix/data")
  indices <- rhdf5::h5read(h5_path, "/matrix/indices")
  indptr <- rhdf5::h5read(h5_path, "/matrix/indptr")
  shape <- rhdf5::h5read(h5_path, "/matrix/shape")
  barcodes <- rhdf5::h5read(h5_path, "/matrix/barcodes")
  genes <- rhdf5::h5read(h5_path, "/matrix/features/name")
  mat <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = as.integer(shape)
  )
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  as(mat, "dgCMatrix")
}

read_tissue_positions <- function(spatial_dir) {
  f <- file.path(spatial_dir, "tissue_positions.csv")
  header <- TRUE
  if (!file.exists(f)) {
    f <- file.path(spatial_dir, "tissue_positions_list.csv")
    header <- FALSE
  }
  if (!file.exists(f)) stop("Could not find tissue positions file in ", spatial_dir)
  pos <- read.csv(f, header = header, stringsAsFactors = FALSE)
  if (!header) {
    colnames(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }
  pos
}

collapse_counts_by_symbol <- function(counts) {
  symbols <- toupper(rownames(counts))
  ok <- !is.na(symbols) & nzchar(symbols)
  counts <- counts[ok, , drop = FALSE]
  symbols <- symbols[ok]
  groups <- factor(symbols, levels = sort(unique(symbols)))
  G <- Matrix::sparse.model.matrix(~ 0 + groups)
  out <- t(G) %*% counts
  rownames(out) <- levels(groups)
  colnames(out) <- colnames(counts)
  as(out, "dgCMatrix")
}

prepare_full_dataset <- function(config, min_detected_fraction = 0.01, scale_factor = 1e4) {
  counts_raw <- read_10x_h5(file.path(config$root, "filtered_feature_bc_matrix.h5"))
  pos <- read_tissue_positions(file.path(config$root, "spatial"))
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(counts_raw), pos$barcode)
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]
  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)
  counts <- collapse_counts_by_symbol(counts_raw)
  min_spots <- max(20L, ceiling(min_detected_fraction * ncol(counts)))
  counts <- counts[Matrix::rowSums(counts > 0) >= min_spots, , drop = FALSE]
  norm <- counts %*% Matrix::Diagonal(x = scale_factor / pmax(lib_size, 1))
  norm@x <- log1p(norm@x)
  Y <- t(as.matrix(norm))
  keep_var <- apply(Y, 2L, stats::var) > 1e-8
  Y <- Y[, keep_var, drop = FALSE]
  coords <- as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  coords_scaled <- scale(coords)
  coords_scaled <- sweep(coords_scaled, 2L, apply(coords_scaled, 2L, min), "-")
  coords_scaled <- sweep(coords_scaled, 2L, apply(coords_scaled, 2L, max), "/")
  list(
    Y = Y,
    coords = coords_scaled,
    spot_metadata = data.frame(
      barcode = common,
      lib_size = as.numeric(lib_size),
      n_detected = as.numeric(n_detected),
      pos,
      row.names = NULL
    )
  )
}

select_reference_markers <- function(ref_log, marker_per_type, min_marker_score) {
  markers <- character()
  rows <- list()
  for (ct in colnames(ref_log)) {
    others <- setdiff(colnames(ref_log), ct)
    score <- ref_log[, ct] - apply(ref_log[, others, drop = FALSE], 1L, max)
    ok <- is.finite(score) & score >= min_marker_score & ref_log[, ct] > 0.2
    ord <- order(score, ref_log[, ct], decreasing = TRUE)
    keep <- rownames(ref_log)[ord[ok[ord]]]
    keep <- head(keep, marker_per_type)
    markers <- unique(c(markers, keep))
    rows[[ct]] <- data.frame(
      cell_type = ct,
      gene = keep,
      marker_score = score[keep],
      ref_logcpm = ref_log[keep, ct],
      stringsAsFactors = FALSE
    )
  }
  list(markers = markers, table = do.call(rbind, rows))
}

project_simplex_rows <- function(W) {
  U <- t(apply(W, 1L, sort, decreasing = TRUE))
  cssv <- t(apply(U, 1L, cumsum)) - 1
  jj <- matrix(seq_len(ncol(W)), nrow = nrow(W), ncol = ncol(W), byrow = TRUE)
  ok <- U - cssv / jj > 0
  rho <- rowSums(ok)
  theta <- cssv[cbind(seq_len(nrow(W)), pmax(rho, 1L))] / pmax(rho, 1L)
  pmax(W - theta, 0)
}

fit_simplex_deconv <- function(Y, S, max_iter = 350L, tol = 1e-7) {
  n <- nrow(Y)
  k <- nrow(S)
  W <- matrix(1 / k, nrow = n, ncol = k)
  colnames(W) <- rownames(S)
  SSt <- tcrossprod(S)
  YSt <- Y %*% t(S)
  eig <- eigen(SSt, symmetric = TRUE, only.values = TRUE)$values
  step <- 1 / max(eig, .Machine$double.eps)
  prev <- Inf
  for (iter in seq_len(max_iter)) {
    grad <- W %*% SSt - YSt
    W <- project_simplex_rows(W - step * grad)
    if (iter %% 25L == 0L || iter == max_iter) {
      resid <- W %*% S - Y
      obj <- sum(resid^2) / length(resid)
      message("    deconv iter ", iter, "/", max_iter, " mse=", signif(obj, 5))
      if (is.finite(prev) && abs(prev - obj) < tol * max(1, prev)) break
      prev <- obj
    }
  }
  W
}

plot_spatial <- function(df, value_col, title, out_file) {
  p <- ggplot(df, aes(x = x, y = y, color = .data[[value_col]])) +
    geom_point(size = 0.55, alpha = 0.9) +
    scale_color_viridis(option = "magma") +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, color = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
  ggsave(out_file, p, width = 5.2, height = 4.6, dpi = 180)
}

make_background_factors <- function(Y, Z, hallmark_genes, n_factors, n_genes) {
  bg_genes <- setdiff(colnames(Y), hallmark_genes)
  Y_bg <- Y[, bg_genes, drop = FALSE]
  vars <- apply(Y_bg, 2L, stats::var)
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(n_genes, length(vars)))]
  Y_bg <- Y_bg[, keep, drop = FALSE]
  X <- cbind(1, scale(Z))
  fit <- qr(X)
  R <- qr.resid(fit, Y_bg)
  R <- scale(R, center = TRUE, scale = TRUE)
  R[!is.finite(R)] <- 0
  s <- svd(R, nu = n_factors, nv = 0)
  V <- s$u[, seq_len(n_factors), drop = FALSE] %*% diag(s$d[seq_len(n_factors)], n_factors)
  colnames(V) <- paste0("bg_factor_", seq_len(n_factors))
  V
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
ref <- readRDS(cfg$reference_rds)
ref_log <- ref$log_cpm
rownames(ref_log) <- toupper(rownames(ref_log))
message("Loaded Wu reference with ", nrow(ref_log), " genes and ", ncol(ref_log), " cell types")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- unique(toupper(hallmark$gene_symbol))
marker_info <- select_reference_markers(ref_log, cfg$marker_per_type, cfg$min_marker_score)
message("Selected ", length(marker_info$markers), " reference marker genes")
write.csv(marker_info$table, file.path(cfg$out_root, "wu_reference_markers.csv"), row.names = FALSE)

for (dataset_id in cfg$datasets) {
  if (!dataset_id %in% names(dataset_configs)) stop("Unknown dataset: ", dataset_id)
  config <- dataset_configs[[dataset_id]]
  message("\n=== Deconvolving ", config$label, " ===")
  dat <- prepare_full_dataset(config)
  message("Prepared full expression matrix: ", nrow(dat$Y), " spots x ", ncol(dat$Y), " genes")
  common_genes <- intersect(intersect(colnames(dat$Y), rownames(ref_log)), marker_info$markers)
  if (length(common_genes) < 50L) stop("Too few marker genes overlap: ", length(common_genes))
  message("Marker overlap: ", length(common_genes), " genes")

  S <- t(ref_log[common_genes, , drop = FALSE])
  Y <- dat$Y[, common_genes, drop = FALSE]
  pooled <- rbind(S, Y)
  center <- colMeans(pooled)
  scalev <- apply(pooled, 2L, stats::sd)
  scalev <- pmax(scalev, sqrt(.Machine$double.eps))
  S_scaled <- sweep(sweep(S, 2L, center, "-"), 2L, scalev, "/")
  Y_scaled <- sweep(sweep(Y, 2L, center, "-"), 2L, scalev, "/")

  message("Fitting simplex deconvolution")
  props <- fit_simplex_deconv(Y_scaled, S_scaled, max_iter = cfg$deconv_iter)
  rownames(props) <- dat$spot_metadata$barcode
  recon <- props %*% S_scaled
  residual_mse <- rowMeans((Y_scaled - recon)^2)
  spot_cor <- vapply(seq_len(nrow(Y_scaled)), function(i) {
    suppressWarnings(stats::cor(Y_scaled[i, ], recon[i, ]))
  }, numeric(1))

  message("Computing background factors from non-Hallmark genes")
  bg <- make_background_factors(
    Y = dat$Y,
    Z = props,
    hallmark_genes = hallmark_genes,
    n_factors = cfg$n_bg_factors,
    n_genes = cfg$n_bg_genes
  )
  rownames(bg) <- dat$spot_metadata$barcode

  out_dir <- file.path(cfg$out_root, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(barcode = rownames(props), props, check.names = FALSE),
            file.path(out_dir, "cell_type_proportions.csv"), row.names = FALSE)
  write.csv(data.frame(barcode = rownames(bg), bg, check.names = FALSE),
            file.path(out_dir, "background_factors.csv"), row.names = FALSE)
  write.csv(data.frame(
    barcode = rownames(props),
    residual_mse = residual_mse,
    reconstruction_cor = spot_cor
  ), file.path(out_dir, "deconvolution_qc.csv"), row.names = FALSE)
  saveRDS(list(
    proportions = props,
    background_factors = bg,
    qc = data.frame(residual_mse = residual_mse, reconstruction_cor = spot_cor),
    marker_genes = common_genes,
    cell_types = colnames(props),
    reference = cfg$reference_rds,
    dataset = dataset_id
  ), file.path(out_dir, "wu_deconvolution_result.rds"))

  md <- dat$spot_metadata
  plot_df <- data.frame(
    barcode = md$barcode,
    x = md$pxl_col_in_fullres,
    y = md$pxl_row_in_fullres,
    props,
    residual_mse = residual_mse,
    reconstruction_cor = spot_cor,
    check.names = FALSE
  )
  for (ct in colnames(props)) {
    plot_spatial(plot_df, ct, paste(dataset_id, ct), file.path(out_dir, paste0("celltype_", ct, ".png")))
  }
  for (bf in colnames(bg)) {
    plot_df[[bf]] <- bg[, bf]
    plot_spatial(plot_df, bf, paste(dataset_id, bf), file.path(out_dir, paste0(bf, ".png")))
  }
  plot_spatial(plot_df, "residual_mse", paste(dataset_id, "deconvolution MSE"), file.path(out_dir, "deconvolution_mse.png"))
  plot_spatial(plot_df, "reconstruction_cor", paste(dataset_id, "deconvolution correlation"), file.path(out_dir, "deconvolution_correlation.png"))

  summary_tbl <- data.frame(
    cell_type = colnames(props),
    mean_proportion = colMeans(props),
    median_proportion = apply(props, 2L, stats::median),
    max_proportion = apply(props, 2L, max),
    stringsAsFactors = FALSE
  )
  write.csv(summary_tbl, file.path(out_dir, "cell_type_summary.csv"), row.names = FALSE)
  writeLines(c(
    paste0("# Wu Deconvolution - ", dataset_id),
    "",
    paste0("- Spots: ", nrow(props)),
    paste0("- Marker genes used: ", length(common_genes)),
    paste0("- Mean reconstruction correlation: ", signif(mean(spot_cor, na.rm = TRUE), 4)),
    paste0("- Median reconstruction correlation: ", signif(stats::median(spot_cor, na.rm = TRUE), 4)),
    "",
    "cell_type | mean | median | max",
    "--- | ---: | ---: | ---:",
    apply(summary_tbl, 1L, function(r) paste(r, collapse = " | "))
  ), file.path(out_dir, "deconvolution_summary.md"))
  message("Wrote deconvolution outputs to ", out_dir)
  print(summary_tbl[order(-summary_tbl$mean_proportion), ], row.names = FALSE)
}
