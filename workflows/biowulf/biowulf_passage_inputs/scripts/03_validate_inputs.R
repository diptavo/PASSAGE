#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "/data/DCEG_Dutta/PASSAGE_production_inputs_20260825"
suppressPackageStartupMessages(library(Matrix))

task_files <- sort(list.files(file.path(root, "metadata", "sample_tasks"), pattern = "[.]csv$", full.names = TRUE))
if (length(task_files) != 4L) stop("Expected four sample task manifests, found ", length(task_files))
manifest <- do.call(rbind, lapply(task_files, read.csv, stringsAsFactors = FALSE))
manifest <- manifest[order(manifest$task_id), , drop = FALSE]

pathway_file <- file.path(root, "pathways", "msigdb_human_pathways_filtered.rds")
pathway_meta_file <- file.path(root, "pathways", "msigdb_human_pathways_filtered_metadata.csv")
pathways <- readRDS(pathway_file)
pathway_meta <- read.csv(pathway_meta_file, stringsAsFactors = FALSE)
if (!is.list(pathways) || !length(pathways)) stop("Pathway file is empty or malformed")

checks <- list()
cell_type_rows <- list()
covariate_rows <- list()
sample_rows <- list()

add_check <- function(sample, check, passed, value, requirement) {
  data.frame(
    sample = sample,
    check = check,
    passed = isTRUE(passed),
    value = as.character(value),
    requirement = requirement,
    stringsAsFactors = FALSE
  )
}

for (ii in seq_len(nrow(manifest))) {
  m <- manifest[ii, , drop = FALSE]
  obj <- readRDS(m$file)
  sample <- obj$sample
  Y <- obj$Y
  counts <- obj$expression$raw_counts
  coords <- obj$coords
  P <- obj$cell_type_proportions
  X <- obj$X
  barcodes <- rownames(Y)
  pathway_overlap <- vapply(pathways, function(g) sum(unique(toupper(g)) %in% colnames(Y)), integer(1))

  sample_checks <- rbind(
    add_check(sample, "schema_version", identical(obj$input_schema_version, "1.0.0"), obj$input_schema_version, "1.0.0"),
    add_check(sample, "expression_dimensions", nrow(Y) == ncol(counts) && ncol(Y) == nrow(counts), paste(dim(Y), collapse = "x"), "spots x genes and aligned raw counts"),
    add_check(sample, "expression_finite", all(is.finite(Y)), sum(!is.finite(Y)), "zero non-finite values"),
    add_check(sample, "expression_barcodes", identical(barcodes, colnames(counts)), sum(barcodes == colnames(counts)), "all barcodes aligned"),
    add_check(sample, "coordinate_dimensions", nrow(coords) == nrow(Y) && ncol(coords) == 2L, paste(dim(coords), collapse = "x"), "n_spots x 2"),
    add_check(sample, "coordinates_finite", all(is.finite(coords)), sum(!is.finite(coords)), "zero non-finite values"),
    add_check(sample, "proportion_dimensions", nrow(P) == nrow(Y) && ncol(P) >= 5L, paste(dim(P), collapse = "x"), "n_spots x at least 5 cell types"),
    add_check(sample, "proportion_barcodes", identical(rownames(P), barcodes), sum(rownames(P) == barcodes), "all barcodes aligned"),
    add_check(sample, "proportion_range", min(P) >= -1e-10 && max(P) <= 1 + 1e-10, paste(signif(range(P), 5), collapse = ","), "within [0,1]"),
    add_check(sample, "proportion_sum", max(abs(rowSums(P) - 1)) < 1e-6, signif(max(abs(rowSums(P) - 1)), 5), "maximum deviation < 1e-6"),
    add_check(sample, "design_dimensions", nrow(X) == nrow(Y), paste(dim(X), collapse = "x"), "n_spots x covariates"),
    add_check(sample, "design_full_rank", qr(X)$rank == ncol(X), paste(qr(X)$rank, ncol(X), sep = "/"), "rank equals number of columns"),
    add_check(sample, "design_intercept", "intercept" %in% colnames(X), paste(colnames(X), collapse = ";"), "contains intercept"),
    add_check(sample, "pathway_coverage", sum(pathway_overlap >= 15L) >= 1000L, sum(pathway_overlap >= 15L), "at least 1000 pathways with >=15 observed genes"),
    add_check(sample, "deconvolution_fit", median(obj$deconvolution_diagnostics$cosine, na.rm = TRUE) >= 0.45, signif(median(obj$deconvolution_diagnostics$cosine, na.rm = TRUE), 5), "median cosine >= 0.45; preferred >= 0.50")
  )
  checks[[ii]] <- sample_checks

  cell_type_rows[[ii]] <- data.frame(
    sample = sample,
    cohort = obj$cohort,
    cell_type = colnames(P),
    mean_fraction = colMeans(P),
    median_fraction = apply(P, 2L, median),
    q95_fraction = apply(P, 2L, quantile, probs = 0.95),
    zero_fraction = colMeans(P < 1e-8),
    composition_reference = colnames(P) == obj$cell_type_reference$composition_reference_omitted_from_design,
    stringsAsFactors = FALSE
  )
  covariate_rows[[ii]] <- data.frame(
    sample = sample,
    covariate = colnames(X),
    class = ifelse(colnames(X) == "intercept", "intercept", ifelse(grepl("^cell_fraction_", colnames(X)), "cell_type_fraction", "technical")),
    stringsAsFactors = FALSE
  )
  sample_rows[[ii]] <- data.frame(
    sample = sample,
    cohort = obj$cohort,
    disease = obj$disease,
    spatial_sample = obj$spatial_sample,
    platform = obj$platform,
    source = obj$source,
    reference_id = obj$cell_type_reference$reference_id,
    reference_url = obj$cell_type_reference$source_url,
    n_spots = nrow(Y),
    n_genes = ncol(Y),
    n_cell_types = ncol(P),
    n_covariates = ncol(X),
    pathway_sets_with_15_genes = sum(pathway_overlap >= 15L),
    bundle_file = m$file,
    md5 = unname(tools::md5sum(m$file)),
    stringsAsFactors = FALSE
  )

  component_dir <- file.path(root, "components", sample)
  dir.create(component_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(cbind(barcode = barcodes, as.data.frame(coords)), gzfile(file.path(component_dir, "spatial_coordinates.csv.gz")), row.names = FALSE)
  write.csv(cbind(barcode = rownames(P), as.data.frame(P, check.names = FALSE)), gzfile(file.path(component_dir, "cell_type_proportions.csv.gz")), row.names = FALSE)
  write.csv(cbind(barcode = rownames(X), as.data.frame(X, check.names = FALSE)), gzfile(file.path(component_dir, "adjustment_covariates.csv.gz")), row.names = FALSE)
  write.csv(obj$gene_data, gzfile(file.path(component_dir, "gene_metadata.csv.gz")), row.names = FALSE)
  write.csv(obj$deconvolution_diagnostics, gzfile(file.path(component_dir, "deconvolution_diagnostics.csv.gz")), row.names = FALSE)
}

checks <- do.call(rbind, checks)
cell_types <- do.call(rbind, cell_type_rows)
covariates <- do.call(rbind, covariate_rows)
samples <- do.call(rbind, sample_rows)
write.csv(manifest, file.path(root, "metadata", "input_bundle_manifest.csv"), row.names = FALSE)
write.csv(samples, file.path(root, "metadata", "sample_metadata.csv"), row.names = FALSE)
write.csv(checks, file.path(root, "metadata", "input_validation.csv"), row.names = FALSE)
write.csv(cell_types, file.path(root, "metadata", "cell_type_proportion_summary.csv"), row.names = FALSE)
write.csv(covariates, file.path(root, "metadata", "adjustment_covariate_manifest.csv"), row.names = FALSE)

collection_col <- if ("collection" %in% names(pathway_meta)) "collection" else names(pathway_meta)[[1L]]
pathway_summary <- as.data.frame(table(collection = pathway_meta[[collection_col]]), stringsAsFactors = FALSE)
names(pathway_summary)[[2L]] <- "n_pathways"
write.csv(pathway_summary, file.path(root, "metadata", "pathway_collection_summary.csv"), row.names = FALSE)

all_passed <- all(checks$passed)
readme <- c(
  "# PASSAGE production inputs: kidney and breast cancer",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Validation: ", if (all_passed) "PASS" else "FAIL"),
  "",
  "## Contents",
  "",
  "1. Expression: raw sparse counts and log1p CPM-10k spot-by-gene matrix inside each RDS bundle.",
  "2. Coordinates: pixel x/y coordinates inside each bundle and components/<sample>/spatial_coordinates.csv.gz.",
  "3. Pathways: filtered human MSigDB collections in pathways/msigdb_human_pathways_filtered.rds.",
  "4. Sample metadata: metadata/sample_metadata.csv and the sample_metadata field in each bundle.",
  "5. Cell types: reference-informed nonnegative fractions inside each bundle and components/<sample>/cell_type_proportions.csv.gz.",
  "6. Covariates: technical QC plus K-1 composition terms inside each bundle and components/<sample>/adjustment_covariates.csv.gz.",
  "",
  "## Primary model",
  "",
  "The primary adjustment design contains an intercept, standardized log library size, detected-gene count, mitochondrial fraction, ribosomal fraction, and all but one cell-type fraction. The omitted composition reference is recorded per sample. This avoids the exact linear dependence created by fractions summing to one.",
  "",
  "Marker-score covariates are retained only as a sensitivity design; they are not labeled as proportions.",
  "Cell-type fractions are relative reference-informed estimates. A median reconstruction cosine of 0.45 is required and 0.50 is preferred; sections below the preferred target should be repeated with the marker-score sensitivity design.",
  "",
  "## Samples",
  "",
  paste0("- ", samples$sample, ": ", samples$n_spots, " spots, ", samples$n_genes, " genes, ", samples$n_cell_types, " cell types, reference ", samples$reference_id),
  "",
  paste0("Validation checks passed: ", sum(checks$passed), "/", nrow(checks))
)
writeLines(readme, file.path(root, "README.md"))
writeLines(if (all_passed) "PASS" else "FAIL", file.path(root, "metadata", "validation_status.txt"))
if (!all_passed) {
  failed <- checks[!checks$passed, , drop = FALSE]
  print(failed)
  stop("Input validation failed for ", nrow(failed), " checks")
}
