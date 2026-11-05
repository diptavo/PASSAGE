#!/usr/bin/env Rscript
# Export one frozen prepared bundle without changing its analysis arrays.
# Rscript export_prepared.R INPUT.rds OUTPUT.rds SOURCE_MANIFEST.csv SAMPLE_ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Supply INPUT.rds OUTPUT.rds SOURCE_MANIFEST.csv SAMPLE_ID")
suppressPackageStartupMessages(library(PASSAGE))
suppressPackageStartupMessages(library(Matrix))
stopifnot(packageVersion("PASSAGE") >= "0.6.1")
sources <- read.csv(args[3L], stringsAsFactors = FALSE, check.names = FALSE)
meta <- sources[sources$sample == args[4L], , drop = FALSE]
stopifnot(nrow(meta) == 1L, !file.exists(args[2L]))
source_md5 <- unname(tools::md5sum(args[1L]))
stopifnot(identical(source_md5, meta$prepared_input_md5))
b <- passage_validate_bundle(args[1L])
stopifnot(identical(b$sample, meta$sample), nrow(b$Y) == meta$n_spots,
          ncol(b$Y) == meta$n_genes, ncol(b$X) == meta$n_covariates,
          length(b$pathway_sets) == 50L)

# Preserve the complete prepared expression, design, counts and annotations.
# Discard duplicate expression arrays and machine-specific provenance fields.
fields <- c("input_schema_version", "sample", "cohort", "donor", "disease",
            "spatial_sample", "platform", "sample_metadata", "Y", "coords", "X",
            "spot_data", "pathway_sets", "cell_type_proportions",
            "deconvolution_diagnostics", "reference_status", "dropped_covariates",
            "pathology_annotations", "pathology_annotation_coverage", "qc",
            "adjustment_covariates", "primary_celltype_covariate_columns",
            "detected_celltype_covariate_columns")
out <- b[intersect(fields, names(b))]
out$raw_counts <- b$expression$raw_counts
if (is.null(out$raw_counts)) stop("Prepared bundle has no raw count matrix")
out$normalization <- "log1p(10000 * counts / original spot library size before gene filtering)"
out$gene_data <- data.frame(gene_symbol = colnames(out$Y), row.names = colnames(out$Y))
if (!is.null(b$gene_data) && "gene_symbol" %in% names(b$gene_data)) {
  ii <- match(colnames(out$Y), b$gene_data$gene_symbol)
  if (anyNA(ii) || anyDuplicated(b$gene_data$gene_symbol)) stop("Ambiguous gene metadata")
  out$gene_data <- b$gene_data[ii, , drop = FALSE]
  rownames(out$gene_data) <- colnames(out$Y)
}
if (!is.null(b$cell_type_reference)) {
  keys <- c("reference_id", "source_url", "method", "fraction_estimation", "genes_used",
            "composition_reference_omitted_from_design", "omitted")
  out$cell_type_reference <- b$cell_type_reference[intersect(keys, names(b$cell_type_reference))]
  if (is.null(out$cell_type_reference$source_url)) {
    out$cell_type_reference$source_url <- if (meta$cohort == "breast") {
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078"
    } else "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630"
  }
}
out$sample_metadata$tissue_context <- meta$tissue_context
out$sample_metadata$source_dataset <- meta$source_dataset
out$source_provenance <- list(
  source_dataset = meta$source_dataset, creators = meta$creators,
  source_url = meta$source_url, source_matrix_url = meta$source_matrix_url,
  license = meta$source_license, license_url = meta$license_url,
  prepared_input_date = "2026-10-22", prepared_input_md5 = source_md5,
  export_date = "2026-11-05", PASSAGE_validation_version = "0.6.1",
  export_changes = paste("All prepared spots, genes, counts, coordinates, design columns,",
    "50 Hallmark sets, composition estimates and existing pathology annotations retained.",
    "Removed duplicate arrays and machine-specific paths; clarified source tissue context.",
    "No additional filtering, normalization, or model fitting during export."),
  gene_set_source = "MSigDB human Hallmark collection",
  gene_set_snapshot = "Frozen 2026-10-22; exact MSigDB release was not recorded",
  gene_set_license = "CC-BY-4.0",
  gene_set_license_url = "https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp",
  gene_set_copyright = paste("Copyright 2004-2025 Broad Institute, Inc.,",
    "Massachusetts Institute of Technology, and Regents of the University of California")
)
out <- passage_validate_bundle(out)
stopifnot(identical(out$Y, b$Y), identical(out$coords, b$coords), identical(out$X, b$X),
          identical(out$raw_counts, b$expression$raw_counts),
          identical(out$pathway_sets, b$pathway_sets),
          identical(out$cell_type_proportions, b$cell_type_proportions),
          identical(out$pathology_annotations, b$pathology_annotations),
          identical(rownames(out$raw_counts), colnames(out$Y)),
          identical(colnames(out$raw_counts), rownames(out$Y)),
          identical(out$spot_data$barcode, rownames(out$Y)),
          qr(out$X)$rank == ncol(out$X))
if (!is.null(out$cell_type_proportions)) {
  stopifnot(identical(rownames(out$cell_type_proportions), rownames(out$Y)),
            all(is.finite(out$cell_type_proportions)),
            min(out$cell_type_proportions) >= -1e-8,
            max(abs(rowSums(out$cell_type_proportions) - 1)) < 1e-6)
}
if (!is.null(out$pathology_annotations)) {
  stopifnot(identical(out$pathology_annotations$barcode, rownames(out$Y)))
}
dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, args[2L], compress = "xz", version = 3)
check <- passage_validate_bundle(args[2L])
stopifnot(identical(check, out), file.info(args[2L])$size < 100 * 1024^2)
cat("Export verified:", meta$sample, nrow(out$Y), "spots;", ncol(out$Y), "genes\n")
print(tools::md5sum(args[2L]))
