# Run SpaPath on two public 10x Visium breast cancer datasets
# with MSigDB Hallmark pathways.
#
# Expected files:
#   data/raw/<dataset>/filtered_feature_bc_matrix.h5
#   data/raw/<dataset>/spatial/tissue_positions_list.csv
#
# Usage:
#   cd /path/to/PASSAGE
#   Rscript scripts/run_10x_breast_hallmark.R

suppressPackageStartupMessages({
  library(Matrix)
  library(msigdbr)
  library(rhdf5)
})

source(file.path("R", "spapath.R"))

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

get_hallmark_pathways <- function() {
  hall <- msigdbr(species = "Homo sapiens", collection = "H")
  hall$gene_symbol <- toupper(hall$gene_symbol)
  lapply(split(hall$gene_symbol, hall$gs_name), unique)
}

read_10x_h5 <- function(h5_path) {
  data <- rhdf5::h5read(h5_path, "/matrix/data")
  indices <- rhdf5::h5read(h5_path, "/matrix/indices")
  indptr <- rhdf5::h5read(h5_path, "/matrix/indptr")
  shape <- rhdf5::h5read(h5_path, "/matrix/shape")
  barcodes <- rhdf5::h5read(h5_path, "/matrix/barcodes")
  genes <- rhdf5::h5read(h5_path, "/matrix/features/name")

  mat <- sparseMatrix(
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
  if (!file.exists(f)) {
    stop("Could not find tissue positions file in ", spatial_dir)
  }
  pos <- read.csv(f, header = header, stringsAsFactors = FALSE)
  if (!header) {
    colnames(pos) <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  }
  pos
}

collapse_counts_by_symbol <- function(counts) {
  symbols <- toupper(rownames(counts))
  ok <- !is.na(symbols) & nzchar(symbols)
  counts <- counts[ok, , drop = FALSE]
  symbols <- symbols[ok]
  groups <- factor(symbols, levels = sort(unique(symbols)))
  G <- sparse.model.matrix(~ 0 + groups)
  out <- t(G) %*% counts
  rownames(out) <- levels(groups)
  colnames(out) <- colnames(counts)
  as(out, "dgCMatrix")
}

prepare_dataset <- function(config,
                            hallmark_pathways,
                            min_detected_fraction = 0.01,
                            scale_factor = 1e4) {
  h5_path <- file.path(config$root, "filtered_feature_bc_matrix.h5")
  spatial_dir <- file.path(config$root, "spatial")

  counts_raw <- read_10x_h5(h5_path)
  pos <- read_tissue_positions(spatial_dir)
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]

  common <- intersect(colnames(counts_raw), pos$barcode)
  if (length(common) == 0) {
    stop("No overlapping barcodes between counts and spatial positions")
  }
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]

  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)

  counts <- collapse_counts_by_symbol(counts_raw)
  hallmark_genes <- sort(unique(unlist(hallmark_pathways, use.names = FALSE)))
  keep_hallmark <- intersect(hallmark_genes, rownames(counts))
  counts <- counts[keep_hallmark, , drop = FALSE]

  min_spots <- max(20, ceiling(min_detected_fraction * ncol(counts)))
  detected <- Matrix::rowSums(counts > 0)
  counts <- counts[detected >= min_spots, , drop = FALSE]

  norm <- counts %*% Diagonal(x = scale_factor / pmax(lib_size, 1))
  norm@x <- log1p(norm@x)
  Y <- t(as.matrix(norm))

  keep_var <- apply(Y, 2, stats::var) > 1e-8
  Y <- Y[, keep_var, drop = FALSE]

  coords <- as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  coords <- scale(coords)
  coords <- sweep(coords, 2, apply(coords, 2, min), "-")
  coords <- sweep(coords, 2, apply(coords, 2, max), "/")

  X <- cbind(
    log_umi = as.numeric(scale(log1p(lib_size))),
    log_detected = as.numeric(scale(log1p(n_detected)))
  )
  rownames(X) <- colnames(counts_raw)

  pathways <- lapply(hallmark_pathways, intersect, y = colnames(Y))
  pathways <- pathways[lengths(pathways) >= 5]

  list(
    Y = Y,
    coords = coords,
    X = X,
    pathways = pathways,
    spot_metadata = data.frame(
      barcode = common,
      lib_size = as.numeric(lib_size),
      n_detected = as.numeric(n_detected),
      pos,
      row.names = NULL
    ),
    n_spots = nrow(Y),
    n_features = ncol(Y),
    n_pathways = length(pathways)
  )
}

flatten_result <- function(res) {
  x <- res$results
  x$driver_genes <- vapply(x$driver_genes, paste, collapse = ";", FUN.VALUE = character(1))
  x$driver_weights <- vapply(x$driver_weights, function(w) {
    paste(signif(w, 4), collapse = ";")
  }, FUN.VALUE = character(1))
  x$driver_scores <- vapply(x$driver_scores, function(z) {
    paste(signif(z, 4), collapse = ";")
  }, FUN.VALUE = character(1))
  x
}

write_outputs <- function(dataset_id, config, dat, fit) {
  out_dir <- file.path("results", dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(fit, file.path(out_dir, "spapath_hallmark_result.rds"))
  write.csv(flatten_result(fit), file.path(out_dir, "spapath_hallmark_pathways.csv"),
            row.names = FALSE)
  write.csv(fit$feature_scores, file.path(out_dir, "spapath_hallmark_feature_scores.csv"),
            row.names = FALSE)
  write.csv(dat$spot_metadata, file.path(out_dir, "spot_metadata.csv"), row.names = FALSE)

  top <- utils::head(flatten_result(fit), 15)
  md <- c(
    paste0("# ", config$label),
    "",
    paste0("- Spots: ", dat$n_spots),
    paste0("- Hallmark features analyzed: ", dat$n_features),
    paste0("- Hallmark pathways tested: ", dat$n_pathways),
    paste0("- Kernel: ", fit$kernel),
    paste0("- Neighbor size: ", fit$m),
    "",
    "## Top Hallmark Pathways",
    "",
    paste(c("pathway", "p_value", "fdr", "eSPVE_any", "q_eff", "best_range", "driver_genes"),
          collapse = " | "),
    paste(rep("---", 7), collapse = " | ")
  )
  rows <- apply(top[, c("pathway", "p_value", "fdr", "eSPVE_any", "q_eff",
                       "best_range", "driver_genes")], 1, function(r) {
    paste(r, collapse = " | ")
  })
  writeLines(c(md, rows), file.path(out_dir, "summary.md"))
}

hallmark_pathways <- get_hallmark_pathways()

for (dataset_id in names(dataset_configs)) {
  config <- dataset_configs[[dataset_id]]
  message("\n=== ", config$label, " ===")
  dat <- prepare_dataset(config, hallmark_pathways)
  message("Prepared ", dat$n_spots, " spots, ", dat$n_features,
          " Hallmark genes, ", dat$n_pathways, " pathways")

  fit <- spapath_test(
    Y = dat$Y,
    coords = dat$coords,
    X = dat$X,
    pathways = dat$pathways,
    m = 20,
    ranges = spapath_default_ranges(dat$coords),
    n_sim = 2000,
    min_pathway_size = 5,
    max_pathway_size = 500,
    seed = 20260523,
    verbose = TRUE
  )

  write_outputs(dataset_id, config, dat, fit)
  print(fit, n = 8)
}
