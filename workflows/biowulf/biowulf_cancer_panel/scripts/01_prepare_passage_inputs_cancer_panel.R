#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)

suppressPackageStartupMessages({
  library(Matrix)
})

out_dir <- file.path(root, "data", "passage_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spatial_manifest <- read.csv(file.path(root, "data", "spatial_manifest.csv"),
                             stringsAsFactors = FALSE)
marker_refs <- read.csv(file.path(root, "refs", "celltype_marker_references.csv"),
                        stringsAsFactors = FALSE)

find_one <- function(root_dir, pattern) {
  hits <- list.files(root_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("__MACOSX|\\.DS_Store", hits)]
  if (!length(hits)) return(NA_character_)
  hits[[1L]]
}

read_table_any <- function(path, header = TRUE) {
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), header = header, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, header = header, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

read_10x <- function(outs_dir) {
  mtx <- find_one(outs_dir, "^matrix[.]mtx([.]gz)?$")
  features <- find_one(outs_dir, "^(features|genes)[.]tsv([.]gz)?$")
  barcodes <- find_one(outs_dir, "^barcodes[.]tsv([.]gz)?$")
  if (any(is.na(c(mtx, features, barcodes)))) {
    stop("Missing 10x matrix/features/barcodes under ", outs_dir)
  }
  M <- Matrix::readMM(if (grepl("[.]gz$", mtx)) gzfile(mtx) else mtx)
  ft <- read_table_any(features, header = FALSE)
  bc <- read_table_any(barcodes, header = FALSE)
  symbols <- if (ncol(ft) >= 2L) ft[[2L]] else ft[[1L]]
  symbols[is.na(symbols) | !nzchar(symbols)] <- ft[[1L]][is.na(symbols) | !nzchar(symbols)]
  rownames(M) <- make.unique(as.character(symbols))
  colnames(M) <- as.character(bc[[1L]])
  list(counts = M, features = ft, barcodes = bc[[1L]])
}

read_positions <- function(outs_dir) {
  pos <- find_one(outs_dir, "^tissue_positions(_list)?[.]csv$")
  if (is.na(pos)) stop("Missing tissue_positions*.csv under ", outs_dir)
  first <- readLines(pos, n = 1L)
  has_header <- grepl("barcode|in_tissue|array_row|pxl", first, ignore.case = TRUE)
  df <- read.csv(pos, header = has_header, stringsAsFactors = FALSE, check.names = FALSE)
  if (!has_header) {
    colnames(df)[seq_len(min(6L, ncol(df)))] <- c("barcode", "in_tissue", "array_row", "array_col",
                                                  "pxl_row_in_fullres", "pxl_col_in_fullres")[seq_len(min(6L, ncol(df)))]
  }
  if (!"barcode" %in% colnames(df)) colnames(df)[[1L]] <- "barcode"
  if (!"in_tissue" %in% colnames(df)) df$in_tissue <- 1L
  x_col <- intersect(c("pxl_col_in_fullres", "imagecol", "array_col"), colnames(df))[1L]
  y_col <- intersect(c("pxl_row_in_fullres", "imagerow", "array_row"), colnames(df))[1L]
  if (is.na(x_col) || is.na(y_col)) stop("Could not identify x/y coordinate columns in ", pos)
  df$x <- as.numeric(df[[x_col]])
  df$y <- as.numeric(df[[y_col]])
  df
}

zscore_rows <- function(M) {
  M <- as.matrix(M)
  mu <- rowMeans(M)
  centered <- M - mu
  v <- rowMeans(centered * centered)
  sd <- sqrt(pmax(v, 1e-8))
  centered / sd
}

marker_scores <- function(log_expr, genes, marker_tbl) {
  gene_upper <- toupper(genes)
  scores <- lapply(seq_len(nrow(marker_tbl)), function(ii) {
    g <- unique(strsplit(marker_tbl$genes[[ii]], ";", fixed = TRUE)[[1L]])
    idx <- which(gene_upper %in% g)
    if (!length(idx)) return(rep(0, ncol(log_expr)))
    Z <- zscore_rows(log_expr[idx, , drop = FALSE])
    as.numeric(colMeans(Z))
  })
  X <- do.call(cbind, scores)
  colnames(X) <- paste0("marker_", make.names(marker_tbl$reference), "_",
                        make.names(marker_tbl$cell_type))
  X[!is.finite(X)] <- 0
  X
}

make_input <- function(row, reference_name) {
  tenx <- read_10x(row$outs_dir)
  pos <- read_positions(row$outs_dir)
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(tenx$counts), pos$barcode)
  if (length(common) < 100L) stop("Too few tissue spots for ", row$spatial_sample, ": ", length(common))
  counts <- tenx$counts[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]
  lib_size <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  keep_gene <- Matrix::rowSums(counts > 0) >= 3L
  counts <- counts[keep_gene, , drop = FALSE]
  genes <- rownames(counts)
  lib_safe <- pmax(lib_size, 1)
  log_expr <- log1p(counts %*% Matrix::Diagonal(x = 1e4 / lib_safe))
  ref <- marker_refs[marker_refs$cancer == row$cancer & marker_refs$reference == reference_name, , drop = FALSE]
  if (!nrow(ref)) stop("No marker reference ", reference_name, " for ", row$cancer)
  X_tech <- data.frame(
    intercept = 1,
    log_library_size = log1p(lib_size[common]),
    detected_genes = detected[common]
  )
  X_cell <- marker_scores(log_expr, genes, ref)
  X <- as.matrix(cbind(X_tech, X_cell))
  sample_id <- paste(row$cancer, row$spatial_sample, reference_name, sep = "__")
  sample_id <- gsub("[^A-Za-z0-9_.-]+", "_", sample_id)
  obj <- list(
    Y = as.matrix(t(log_expr)),
    coords = as.matrix(pos[, c("x", "y"), drop = FALSE]),
    X = X,
    spot_data = pos,
    gene_data = data.frame(gene_symbol = genes, stringsAsFactors = FALSE),
    assay_name = "log1p_cpm_10k",
    sample = sample_id,
    cancer = row$cancer,
    spatial_sample = row$spatial_sample,
    reference = reference_name,
    source = row$source,
    primary_celltype_covariate_columns = colnames(X_cell),
    detected_celltype_covariate_columns = colnames(X_cell)
  )
  out <- file.path(out_dir, paste0("passage_input_", sample_id, ".rds"))
  saveRDS(obj, out, compress = "gzip")
  data.frame(
    sample = sample_id,
    cancer = row$cancer,
    spatial_sample = row$spatial_sample,
    reference = reference_name,
    n_spots = nrow(obj$Y),
    n_genes = ncol(obj$Y),
    n_covariates = ncol(X),
    has_celltype_covariates = TRUE,
    file = out,
    stringsAsFactors = FALSE
  )
}

ok <- spatial_manifest$matrix_ok & spatial_manifest$spatial_ok &
  spatial_manifest$matrix_untar_ok & spatial_manifest$spatial_untar_ok
spatial_manifest <- spatial_manifest[ok, , drop = FALSE]
if (!nrow(spatial_manifest)) stop("No completed spatial downloads in spatial_manifest.csv")

refs <- c("tumor_tme_broad", "normal_tissue_broad")
manifest <- list()
for (ii in seq_len(nrow(spatial_manifest))) {
  for (rr in refs) {
    message("Preparing ", spatial_manifest$cancer[[ii]], " / ", spatial_manifest$spatial_sample[[ii]],
            " with reference ", rr)
    manifest[[length(manifest) + 1L]] <- make_input(spatial_manifest[ii, ], rr)
  }
}
manifest <- do.call(rbind, manifest)
write.csv(manifest, file.path(out_dir, "passage_input_manifest.csv"), row.names = FALSE)
message("Prepared ", nrow(manifest), " PASSAGE cancer-panel input files in ", out_dir)
