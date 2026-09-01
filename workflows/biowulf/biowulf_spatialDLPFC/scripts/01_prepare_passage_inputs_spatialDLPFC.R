#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
max_samples <- if (length(args) >= 2L) as.integer(args[[2L]]) else NA_integer_

data_dir <- file.path(root, "data", "spatialDLPFC")
out_dir <- file.path(root, "data", "passage_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

spe_path <- file.path(data_dir, "spatialDLPFC_Visium_spe.rds")
if (!file.exists(spe_path)) stop("Missing ", spe_path, ". Run 00_download first.")
spe <- readRDS(spe_path)

assay_name <- if ("logcounts" %in% assayNames(spe)) "logcounts" else "counts"
expr_assay <- assay(spe, assay_name)
count_assay <- if ("counts" %in% assayNames(spe)) assay(spe, "counts") else NULL

cd <- as.data.frame(colData(spe))
sample_col <- intersect(c("sample_id", "sample_name", "sample", "array"),
                        colnames(cd))[1L]
if (is.na(sample_col)) {
  sample_col <- "all_spots"
  cd[[sample_col]] <- "all_spots"
}

coords <- as.data.frame(spatialCoords(spe))
if (ncol(coords) < 2L) {
  coord_cols <- intersect(c("array_col", "array_row", "pxl_col_in_fullres",
                            "pxl_row_in_fullres"), colnames(cd))
  coords <- cd[, coord_cols[seq_len(min(2L, length(coord_cols)))], drop = FALSE]
}
if (ncol(coords) < 2L) stop("Could not find at least two spatial coordinate columns.")
coords <- coords[, seq_len(2L), drop = FALSE]
colnames(coords) <- c("x", "y")

numeric_cd <- vapply(cd, is.numeric, logical(1))
candidate_celltype_cols <- grep(
  "(^prop_|_prop$|proportion|decon|cell2location|rctd|spotlight|card)",
  colnames(cd), ignore.case = TRUE, value = TRUE
)
candidate_celltype_cols <- candidate_celltype_cols[numeric_cd[candidate_celltype_cols]]
if (length(candidate_celltype_cols)) {
  prop_mat <- as.matrix(cd[, candidate_celltype_cols, drop = FALSE])
  ok_range <- colMeans(prop_mat >= -0.05 & prop_mat <= 1.05, na.rm = TRUE) > 0.95
  candidate_celltype_cols <- candidate_celltype_cols[ok_range]
}
write.csv(data.frame(celltype_covariate_column = candidate_celltype_cols),
          file.path(out_dir, "detected_celltype_covariate_columns.csv"),
          row.names = FALSE)
primary_celltype_cols <- grep("^broad_cell2location_", candidate_celltype_cols,
                              value = TRUE)
if (!length(primary_celltype_cols)) {
  primary_celltype_cols <- grep("^broad_spotlight_", candidate_celltype_cols,
                                value = TRUE)
}
if (!length(primary_celltype_cols)) primary_celltype_cols <- candidate_celltype_cols
write.csv(data.frame(celltype_covariate_column = primary_celltype_cols),
          file.path(out_dir, "primary_celltype_covariate_columns.csv"),
          row.names = FALSE)

samples <- unique(cd[[sample_col]])
samples <- samples[!is.na(samples)]
if (is.finite(max_samples)) samples <- head(samples, max_samples)

manifest <- list()
for (s in samples) {
  keep <- cd[[sample_col]] == s
  if (sum(keep) < 100L) next
  expr_sample <- as.matrix(expr_assay[, keep, drop = FALSE])
  if (assay_name == "counts") expr_sample <- log1p(expr_sample)
  if (!is.null(count_assay)) {
    count_sample <- as.matrix(count_assay[, keep, drop = FALSE])
    lib_size <- colSums(count_sample)
    detected <- colSums(count_sample > 0)
  } else {
    lib_size <- colSums(expm1(expr_sample))
    detected <- colSums(expr_sample > 0)
  }
  X_tech <- data.frame(
    intercept = 1,
    log_library_size = log1p(lib_size),
    detected_genes = detected
  )
  X_cell <- NULL
  if (length(primary_celltype_cols)) {
    P <- as.matrix(cd[keep, primary_celltype_cols, drop = FALSE])
    P[!is.finite(P)] <- 0
    P <- pmax(P, 0)
    X_cell <- P[, seq_len(max(1L, ncol(P) - 1L)), drop = FALSE]
    colnames(X_cell) <- paste0("celltype_", make.names(colnames(X_cell)))
  }
  X <- as.matrix(cbind(X_tech, X_cell))
  obj <- list(
    Y = t(expr_sample),
    coords = as.matrix(coords[keep, , drop = FALSE]),
    X = X,
    spot_data = cd[keep, , drop = FALSE],
    gene_data = as.data.frame(rowData(spe)),
    assay_name = assay_name,
    sample = s,
    sample_col = sample_col,
    primary_celltype_covariate_columns = primary_celltype_cols,
    detected_celltype_covariate_columns = candidate_celltype_cols
  )
  nm <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(s))
  out <- file.path(out_dir, paste0("passage_input_", nm, ".rds"))
  saveRDS(obj, out, compress = "gzip")
  manifest[[length(manifest) + 1L]] <- data.frame(
    sample = as.character(s),
    n_spots = sum(keep),
    n_genes = nrow(expr_assay),
    n_covariates = ncol(X),
    has_celltype_covariates = length(primary_celltype_cols) > 0L,
    file = out,
    stringsAsFactors = FALSE
  )
}

manifest <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
write.csv(manifest, file.path(out_dir, "passage_input_manifest.csv"),
          row.names = FALSE)
message("Prepared ", nrow(manifest), " PASSAGE input files in ", out_dir)
