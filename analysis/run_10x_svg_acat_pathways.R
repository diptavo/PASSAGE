# Residualized gene-level SVG kernel tests aggregated to pathways by ACAT.
#
# This is a factor-model-free comparator for PASSAGE. Each gene is tested for
# positive spatial covariance with sparse Matérn kernels over a range grid after
# residualizing the same technical/spatial covariates used in PASSAGE H1.

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
})

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

parse_args <- function(args) {
  cfg <- list(
    out_root = "results/svg_acat_10x",
    pathway_gmt = file.path("msigdb", "h.all.v2023.2.Hs.symbols.gmt"),
    pathway_prefix = "HALLMARK_",
    pathway_label = "HALLMARK",
    datasets = c("Visium_FFPE_Human_Breast_Cancer", "V1_Breast_Cancer_Block_A_Section_1"),
    m = 20L,
    kernels = c("matern12", "matern32", "matern52"),
    n_range = 5L,
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    min_spot_umi = 500,
    min_spot_genes = 200,
    max_mito_fraction = 0.25,
    max_library_quantile = 0.995
  )
  for (arg in args) {
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--pathway-gmt=", arg)) cfg$pathway_gmt <- sub("^--pathway-gmt=", "", arg)
    if (grepl("^--pathway-prefix=", arg)) cfg$pathway_prefix <- sub("^--pathway-prefix=", "", arg)
    if (grepl("^--pathway-label=", arg)) cfg$pathway_label <- sub("^--pathway-label=", "", arg)
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--kernels=", arg)) cfg$kernels <- strsplit(sub("^--kernels=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--n-range=", arg)) cfg$n_range <- as.integer(sub("^--n-range=", "", arg))
    if (grepl("^--min-pathway-size=", arg)) cfg$min_pathway_size <- as.integer(sub("^--min-pathway-size=", "", arg))
    if (grepl("^--max-pathway-size=", arg)) cfg$max_pathway_size <- as.integer(sub("^--max-pathway-size=", "", arg))
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

read_gmt_pathways <- function(gmt_path, prefix = NULL) {
  if (!file.exists(gmt_path)) stop("GMT file not found: ", gmt_path)
  lines <- readLines(gmt_path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) next
    nm <- fields[[1L]]
    if (!is.null(prefix) && nzchar(prefix) && !startsWith(nm, prefix)) next
    out[[nm]] <- toupper(unique(fields[-c(1L, 2L)]))
  }
  out
}

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

prepare_dataset <- function(config, pathways, cfg, scale_factor = 1e4) {
  counts_raw <- read_10x_h5(file.path(config$root, "filtered_feature_bc_matrix.h5"))
  pos <- read_tissue_positions(file.path(config$root, "spatial"))
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(counts_raw), pos$barcode)
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]

  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)
  mito_rows <- grepl("^MT-", toupper(rownames(counts_raw)))
  mito_counts <- if (any(mito_rows)) Matrix::colSums(counts_raw[mito_rows, , drop = FALSE]) else rep(0, ncol(counts_raw))
  mito_fraction <- as.numeric(mito_counts / pmax(lib_size, 1))
  max_lib <- as.numeric(stats::quantile(lib_size, probs = cfg$max_library_quantile, na.rm = TRUE, names = FALSE))
  keep <- lib_size >= cfg$min_spot_umi &
    n_detected >= cfg$min_spot_genes &
    mito_fraction <= cfg$max_mito_fraction &
    lib_size <= max_lib

  counts_raw <- counts_raw[, keep, drop = FALSE]
  pos <- pos[keep, , drop = FALSE]
  lib_size <- lib_size[keep]
  n_detected <- n_detected[keep]
  common <- common[keep]

  counts <- collapse_counts_by_symbol(counts_raw)
  pathway_genes <- sort(unique(unlist(pathways, use.names = FALSE)))
  counts <- counts[intersect(pathway_genes, rownames(counts)), , drop = FALSE]
  min_spots <- max(20L, ceiling(0.01 * ncol(counts)))
  counts <- counts[Matrix::rowSums(counts > 0) >= min_spots, , drop = FALSE]
  norm <- counts %*% Matrix::Diagonal(x = scale_factor / pmax(lib_size, 1))
  norm@x <- log1p(norm@x)
  Y <- t(as.matrix(norm))
  keep_var <- apply(Y, 2L, stats::var) > 1e-8
  Y <- Y[, keep_var, drop = FALSE]

  coords <- as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  coords <- scale(coords)
  coords <- sweep(coords, 2L, apply(coords, 2L, min), "-")
  coords <- sweep(coords, 2L, apply(coords, 2L, max), "/")
  sx <- as.numeric(scale(coords[, 1L]))
  sy <- as.numeric(scale(coords[, 2L]))
  X <- cbind(
    log_umi = as.numeric(scale(log1p(lib_size))),
    log_detected = as.numeric(scale(log1p(n_detected))),
    spatial_x = sx,
    spatial_y = sy,
    spatial_x2 = as.numeric(scale(sx^2)),
    spatial_y2 = as.numeric(scale(sy^2)),
    spatial_xy = as.numeric(scale(sx * sy))
  )
  rownames(X) <- common
  list(
    Y = Y,
    coords = coords,
    X = X,
    spot_metadata = data.frame(barcode = common, lib_size = as.numeric(lib_size), n_detected = as.numeric(n_detected), pos),
    pathways = lapply(pathways, intersect, y = colnames(Y))
  )
}

gene_svg_kernel <- function(Y, coords, X, cfg) {
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  qx <- qr(X)
  R <- passage_residualize_with_qr(Y, qx)
  df <- nrow(R) - qx$rank
  sigma2 <- pmax(colSums(R^2) / max(df, 1L), .Machine$double.eps)
  range_grid <- passage_default_range_grid(coords, n_grid = cfg$n_range, min_frac = 0.04, max_frac = 0.50)
  p_by_kernel <- list()
  q_by_kernel <- list()
  for (kernel in cfg$kernels) {
    p_mat <- matrix(NA_real_, nrow = ncol(Y), ncol = length(range_grid), dimnames = list(colnames(Y), paste0("range", seq_along(range_grid))))
    q_mat <- p_mat
    for (rr in seq_along(range_grid)) {
      vc <- passage_vecchia_precision(coords, range = range_grid[rr], m = cfg$m, kernel = kernel, ordering = "coordinate")
      Ksp <- passage_sparse_covariance_from_vecchia(vc)
      X_ord <- X[vc$ord, , drop = FALSE]
      R_ord <- R[vc$ord, , drop = FALSE]
      moments <- passage_kernel_moments(Ksp, X_ord)
      KR <- Ksp %*% R_ord
      q <- colSums(R_ord * as.matrix(KR))
      mean <- sigma2 * moments$trace_mk
      var <- 2 * sigma2^2 * moments$trace_mkmk
      p <- mapply(passage_satterthwaite_p, q = q, mean = mean, var = var)
      p_mat[, rr] <- p
      q_mat[, rr] <- q
    }
    p_by_kernel[[kernel]] <- p_mat
    q_by_kernel[[kernel]] <- q_mat
  }
  gene_p <- vapply(seq_len(ncol(Y)), function(j) {
    p <- unlist(lapply(p_by_kernel, function(m) m[j, ]), use.names = FALSE)
    passage_acat(p)
  }, numeric(1))
  data.frame(
    gene = colnames(Y),
    svg_acat_p = gene_p,
    svg_acat_fdr = stats::p.adjust(gene_p, method = "BH"),
    stringsAsFactors = FALSE
  )
}

pathway_acat <- function(pathways, gene_svg, min_size, max_size) {
  rows <- lapply(names(pathways), function(nm) {
    genes <- intersect(pathways[[nm]], gene_svg$gene)
    if (length(genes) < min_size || length(genes) > max_size) return(NULL)
    p <- gene_svg$svg_acat_p[match(genes, gene_svg$gene)]
    data.frame(
      pathway = nm,
      pathway_size = length(genes),
      svg_acat_p = passage_acat(p),
      min_gene_p = min(p, na.rm = TRUE),
      median_gene_p = median(p, na.rm = TRUE),
      top_svg_genes = paste(head(genes[order(p)], 25), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  tbl$svg_acat_fdr <- stats::p.adjust(tbl$svg_acat_p, method = "BH")
  tbl[order(tbl$svg_acat_p), , drop = FALSE]
}

run_dataset <- function(dataset_id, pathways, cfg) {
  config <- dataset_configs[[dataset_id]]
  message("\n=== SVG ACAT ", cfg$pathway_label, ": ", config$label, " ===")
  dat <- prepare_dataset(config, pathways, cfg)
  message("Prepared ", nrow(dat$Y), " spots x ", ncol(dat$Y), " genes")
  gene_svg <- gene_svg_kernel(dat$Y, dat$coords, dat$X, cfg)
  path_tbl <- pathway_acat(dat$pathways, gene_svg, cfg$min_pathway_size, cfg$max_pathway_size)
  out_dir <- file.path(cfg$out_root, cfg$pathway_label, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(gene_svg, file.path(out_dir, "gene_svg_kernel_acat.csv"), row.names = FALSE)
  write.csv(path_tbl, file.path(out_dir, "pathway_svg_acat.csv"), row.names = FALSE)
  write.csv(dat$spot_metadata, file.path(out_dir, "spot_metadata.csv"), row.names = FALSE)
  message("Wrote SVG ACAT outputs to ", out_dir)
  print(head(path_tbl, 12), row.names = FALSE)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
pathways <- read_gmt_pathways(cfg$pathway_gmt, cfg$pathway_prefix)
if (length(pathways) == 0L) stop("No pathways loaded")
bad <- setdiff(cfg$datasets, names(dataset_configs))
if (length(bad)) stop("Unknown dataset id(s): ", paste(bad, collapse = ", "))
for (dataset_id in cfg$datasets) run_dataset(dataset_id, pathways, cfg)
