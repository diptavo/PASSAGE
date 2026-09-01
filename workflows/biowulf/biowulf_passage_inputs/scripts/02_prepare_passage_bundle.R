#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "/data/DCEG_Dutta/PASSAGE_production_inputs_20260825"
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
source(file.path(root, "scripts", "reference_common.R"))
suppressPackageStartupMessages(library(nnls))
suppressPackageStartupMessages(library(quadprog))

kidney_root <- "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803"
breast_root <- "/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803"
tasks <- data.frame(
  task_id = 1:4,
  cohort = c("kidney", "kidney", "kidney", "breast"),
  spatial_sample = c("KC1", "KC2", "KC3", "Visium_FFPE_Human_Breast_Cancer"),
  matrix_dir = c(
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC1/filtered_feature_bc_matrix"),
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC2/filtered_feature_bc_matrix"),
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC3/filtered_feature_bc_matrix"),
    file.path(breast_root, "data/spatial/breast/Visium_FFPE_Human_Breast_Cancer/outs")
  ),
  spatial_dir = c(
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC1/spatial"),
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC2/spatial"),
    file.path(kidney_root, "data/spatial/TLS_VISIUM_USZ/TLS_VISIUM_USZ/10x_Visium/KC3/spatial"),
    file.path(breast_root, "data/spatial/breast/Visium_FFPE_Human_Breast_Cancer/outs")
  ),
  disease = c("renal_cell_carcinoma", "renal_cell_carcinoma", "renal_cell_carcinoma", "breast_carcinoma"),
  platform = c("10x_Visium", "10x_Visium", "10x_Visium", "10x_Visium_FFPE"),
  source = c(
    rep("Zenodo record 14620362 TLS Visium kidney cancer dataset", 3L),
    "10x Genomics public Visium FFPE Human Breast Cancer"
  ),
  stringsAsFactors = FALSE
)
if (!task_id %in% tasks$task_id) stop("Invalid task_id: ", task_id)
row <- tasks[tasks$task_id == task_id, , drop = FALSE]

find_one <- function(directory, pattern) {
  hits <- list.files(directory, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("__MACOSX|[.]DS_Store|[.]ipynb_checkpoints", hits)]
  if (!length(hits)) return(NA_character_)
  hits[[1L]]
}

read_10x <- function(directory) {
  mtx <- find_one(directory, "^matrix[.]mtx([.]gz)?$")
  features <- find_one(directory, "^(features|genes)[.]tsv([.]gz)?$")
  barcodes <- find_one(directory, "^barcodes[.]tsv([.]gz)?$")
  if (any(is.na(c(mtx, features, barcodes)))) stop("Incomplete 10x matrix under ", directory)
  M <- readMM(if (grepl("[.]gz$", mtx)) gzfile(mtx) else mtx)
  ft <- read_table_any(features, header = FALSE)
  bc <- read_table_any(barcodes, header = FALSE)[[1L]]
  symbols <- if (ncol(ft) >= 2L) ft[[2L]] else ft[[1L]]
  symbols[is.na(symbols) | !nzchar(symbols)] <- ft[[1L]][is.na(symbols) | !nzchar(symbols)]
  rownames(M) <- clean_symbol(symbols)
  colnames(M) <- as.character(bc)
  list(counts = collapse_sparse_rows(M, rownames(M)), features = ft)
}

read_positions <- function(directory) {
  path <- find_one(directory, "^tissue_positions(_list)?[.]csv$")
  if (is.na(path)) stop("Missing tissue positions under ", directory)
  first <- readLines(path, n = 1L)
  has_header <- grepl("barcode|in_tissue|array_row|pxl", first, ignore.case = TRUE)
  pos <- read.csv(path, header = has_header, stringsAsFactors = FALSE, check.names = FALSE)
  if (!has_header) {
    names(pos)[seq_len(min(6L, ncol(pos)))] <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")[seq_len(min(6L, ncol(pos)))]
  }
  if (!"barcode" %in% names(pos)) names(pos)[[1L]] <- "barcode"
  if (!"in_tissue" %in% names(pos)) pos$in_tissue <- 1L
  x_col <- intersect(c("pxl_col_in_fullres", "imagecol", "array_col"), names(pos))[[1L]]
  y_col <- intersect(c("pxl_row_in_fullres", "imagerow", "array_row"), names(pos))[[1L]]
  pos$x <- as.numeric(pos[[x_col]])
  pos$y <- as.numeric(pos[[y_col]])
  pos[pos$in_tissue == 1L, , drop = FALSE]
}

estimate_fractions <- function(counts, signature, workers = 1L) {
  common <- intersect(rownames(signature), rownames(counts))
  if (length(common) < 25L) stop("Only ", length(common), " reference genes overlap spatial data")
  S <- pmax(as.matrix(signature[common, , drop = FALSE]), 0)
  B <- as.matrix(counts[common, , drop = FALSE])
  B <- sweep(B, 2L, pmax(colSums(B), 1), "/")
  S <- sweep(S, 2L, pmax(colSums(S), 1e-12), "/")
  background <- rowMeans(S)
  background <- background / pmax(sum(background), 1e-12)
  S <- cbind(S, unmodeled_background = background)
  specificity <- apply(S, 1L, max) / pmax(rowMeans(S), 1e-8)
  weight <- sqrt(pmin(pmax(specificity, 1), 10))
  A <- S * weight
  B <- B * weight

  fit_one <- function(ii) {
    b <- B[, ii]
    if (!sum(b) > 0) return(c(rep(1 / ncol(A), ncol(A)), cosine = NA_real_, rmse = NA_real_))
    k <- ncol(A)
    gram <- crossprod(A)
    lambda <- 1e-4 * mean(diag(gram))
    prior <- c(rep(0.8 / (k - 1L), k - 1L), 0.2)
    Dmat <- gram + diag(lambda + 1e-10, k)
    dvec <- as.numeric(crossprod(A, b)) + lambda * prior
    Amat <- cbind(rep(1, k), diag(k))
    fit <- try(quadprog::solve.QP(Dmat, dvec, Amat, c(1, rep(0, k)), meq = 1L), silent = TRUE)
    p <- if (inherits(fit, "try-error")) prior else pmax(fit$solution, 0)
    p <- p / pmax(sum(p), 1e-12)
    fitted <- as.numeric(A %*% p)
    cosine <- sum(fitted * b) / sqrt(pmax(sum(fitted^2) * sum(b^2), 1e-20))
    c(p, cosine = cosine, rmse = sqrt(mean((fitted - b)^2)))
  }
  ans <- if (.Platform$OS.type == "unix" && workers > 1L) {
    parallel::mclapply(seq_len(ncol(B)), fit_one, mc.cores = workers, mc.preschedule = TRUE)
  } else {
    lapply(seq_len(ncol(B)), fit_one)
  }
  ans <- do.call(rbind, ans)
  proportions <- ans[, seq_len(ncol(A)), drop = FALSE]
  colnames(proportions) <- colnames(A)
  rownames(proportions) <- colnames(counts)
  diagnostics <- data.frame(
    barcode = colnames(counts),
    cosine = ans[, "cosine"],
    rmse = ans[, "rmse"],
    max_fraction = apply(proportions, 1L, max),
    effective_cell_types = 1 / rowSums(proportions^2),
    stringsAsFactors = FALSE
  )
  list(proportions = proportions, diagnostics = diagnostics, genes = common, weights = weight)
}

marker_score_matrix <- function(Y, markers) {
  result <- lapply(names(markers), function(cell_type) {
    idx <- which(colnames(Y) %in% markers[[cell_type]])
    if (!length(idx)) return(rep(0, nrow(Y)))
    M <- Y[, idx, drop = FALSE]
    M <- scale(M)
    M[!is.finite(M)] <- 0
    rowMeans(M)
  })
  out <- do.call(cbind, result)
  colnames(out) <- paste0("marker_score_", names(markers))
  out
}

spatial <- read_10x(row$matrix_dir)
pos <- read_positions(row$spatial_dir)
common_spots <- intersect(colnames(spatial$counts), pos$barcode)
if (length(common_spots) < 100L) stop("Too few matched tissue spots: ", length(common_spots))
counts <- spatial$counts[, common_spots, drop = FALSE]
pos <- pos[match(common_spots, pos$barcode), , drop = FALSE]
rownames(pos) <- common_spots

lib_size <- Matrix::colSums(counts)
detected_genes <- Matrix::colSums(counts > 0)
pct_mito <- 100 * Matrix::colSums(counts[grepl("^MT-", rownames(counts)), , drop = FALSE]) / pmax(lib_size, 1)
pct_ribo <- 100 * Matrix::colSums(counts[grepl("^RP[SL][0-9]", rownames(counts)), , drop = FALSE]) / pmax(lib_size, 1)
keep_gene <- Matrix::rowSums(counts > 0) >= 3L
counts <- counts[keep_gene, , drop = FALSE]
log_expr <- log1p(counts %*% Diagonal(x = 1e4 / pmax(lib_size, 1)))
rownames(log_expr) <- rownames(counts)
colnames(log_expr) <- colnames(counts)
Y <- as.matrix(t(log_expr))

reference_file <- file.path(root, "refs", "derived", paste0(row$cohort, if (row$cohort == "breast") "_GSE176078_signature.rds" else "_GSE224630_signature.rds"))
reference <- readRDS(reference_file)
workers <- max(1L, min(8L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))
deconv <- estimate_fractions(counts, reference$signature, workers)
P <- deconv$proportions

spot_data <- pos
spot_data$library_size <- as.numeric(lib_size)
spot_data$detected_genes <- as.numeric(detected_genes)
spot_data$pct_mito <- as.numeric(pct_mito)
spot_data$pct_ribo <- as.numeric(pct_ribo)
spot_data <- cbind(spot_data, as.data.frame(P, check.names = FALSE))

technical_raw <- data.frame(
  log_library_size = log1p(lib_size),
  detected_genes = detected_genes,
  pct_mito = pct_mito,
  pct_ribo = pct_ribo,
  row.names = common_spots,
  check.names = FALSE
)
technical_scaled <- scale(technical_raw)
technical_scaled[!is.finite(technical_scaled)] <- 0
composition_reference <- names(which.max(colMeans(P)))
P_design <- P[, setdiff(colnames(P), composition_reference), drop = FALSE]
P_scaled <- scale(P_design)
P_scaled[!is.finite(P_scaled)] <- 0
colnames(P_scaled) <- paste0("cell_fraction_", colnames(P_scaled))
X_candidate <- cbind(intercept = 1, technical_scaled, P_scaled)
qr_x <- qr(X_candidate)
keep_x <- sort(qr_x$pivot[seq_len(qr_x$rank)])
X <- X_candidate[, keep_x, drop = FALSE]
rownames(X) <- common_spots

marker_scores <- marker_score_matrix(Y, canonical_markers[[row$cohort]])
rownames(marker_scores) <- common_spots
gene_data <- data.frame(gene_symbol = rownames(counts), stringsAsFactors = FALSE)
rownames(gene_data) <- rownames(counts)

sample_id <- paste(row$cohort, row$spatial_sample, sep = "__")
object <- list(
  input_schema_version = "1.0.0",
  sample = sample_id,
  cohort = row$cohort,
  disease = row$disease,
  spatial_sample = row$spatial_sample,
  platform = row$platform,
  source = row$source,
  expression = list(
    raw_counts = counts,
    normalized_spot_by_gene = Y,
    assay_name = "log1p_cpm_10k",
    normalization = "log1p(10000 * counts / spot library size)"
  ),
  Y = Y,
  coordinates = as.matrix(pos[, c("x", "y"), drop = FALSE]),
  coords = as.matrix(pos[, c("x", "y"), drop = FALSE]),
  spot_data = spot_data,
  gene_data = gene_data,
  pathways = list(
    file = file.path(root, "pathways", "msigdb_human_pathways_filtered.rds"),
    metadata_file = file.path(root, "pathways", "msigdb_human_pathways_filtered_metadata.csv")
  ),
  sample_metadata = as.list(row[1L, c("cohort", "spatial_sample", "disease", "platform", "source")]),
  cell_type_proportions = P,
  cell_type_reference = list(
    reference_id = reference$reference_id,
    source_url = reference$source_url,
    method = reference$method,
    fraction_estimation = "simplex-constrained quadratic programming with nonnegative fractions, weak ridge stabilization, and an unmodeled-background component",
    signature_file = reference_file,
    genes_used = deconv$genes,
    composition_reference_omitted_from_design = composition_reference
  ),
  deconvolution_diagnostics = deconv$diagnostics,
  adjustment_covariates = list(
    primary = X,
    technical_raw = technical_raw,
    technical_scaled = technical_scaled,
    composition_full = P,
    composition_design = P_scaled,
    marker_score_sensitivity = marker_scores
  ),
  X = X,
  primary_celltype_covariate_columns = colnames(P_scaled),
  detected_celltype_covariate_columns = colnames(P_scaled)
)

out_dir <- file.path(root, "inputs", row$cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(sample_id, ".rds"))
saveRDS(object, out_file, compress = "gzip")

manifest_row <- data.frame(
  task_id = task_id,
  sample = sample_id,
  cohort = row$cohort,
  disease = row$disease,
  spatial_sample = row$spatial_sample,
  platform = row$platform,
  reference_id = reference$reference_id,
  n_spots = nrow(Y),
  n_genes = ncol(Y),
  n_cell_types = ncol(P),
  n_adjustment_covariates = ncol(X),
  composition_reference = composition_reference,
  median_deconvolution_cosine = median(deconv$diagnostics$cosine, na.rm = TRUE),
  file = out_file,
  stringsAsFactors = FALSE
)
dir.create(file.path(root, "metadata", "sample_tasks"), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest_row, file.path(root, "metadata", "sample_tasks", sprintf("sample_task_%02d.csv", task_id)), row.names = FALSE)
