#!/usr/bin/env Rscript
# Repackage the frozen prepared KC2 input without changing its analysis arrays.
# Usage: Rscript inst/examples/prepare_kidney.R prepared_KC2.rds public_KC2.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_kidney.R INPUT.rds OUTPUT.rds")
if (file.exists(args[2L])) stop("Output already exists: ", args[2L])
suppressPackageStartupMessages(library(PASSAGE))
source_md5 <- unname(tools::md5sum(args[1L]))
if (!identical(source_md5, "ba665dc5727ba41edef9c0dd97658955")) {
  stop("Input is not the frozen KC2 bundle used for this example")
}
b <- passage_validate_bundle(args[1L])
stopifnot(identical(dim(b$Y), c(2166L, 17630L)), ncol(b$X) == 13L,
          identical(b$sample, "kidney__KC2"))
pathways <- c("HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS",
              "HALLMARK_INFLAMMATORY_RESPONSE")
stopifnot(all(pathways %in% names(b$pathway_sets)))
reference <- b$cell_type_reference[c("reference_id", "source_url", "method",
  "fraction_estimation", "genes_used", "composition_reference_omitted_from_design")]
example <- list(
  input_schema_version = "1.1.0", sample = b$sample,
  sample_metadata = b$sample_metadata,
  Y = b$Y, coords = b$coords, X = b$X,
  raw_counts = b$expression$raw_counts,
  normalization = b$expression$normalization,
  pathway_sets = b$pathway_sets[pathways],
  cell_type_proportions = b$cell_type_proportions,
  cell_type_reference = reference,
  spot_qc = b$spot_data[, c("barcode", "library_size", "detected_genes",
                           "pct_mito", "pct_ribo")],
  source_provenance = list(
    dataset = "10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung (5) Cancer with Tertiary Lymphoid Structures",
    creators = c("Sebastian Dawo", "Kalin Nonchev", "Karina Silina"),
    source_url = "https://zenodo.org/records/14620362",
    doi = "10.5281/zenodo.14620362", source_version = "v1, 2025-01-09",
    license = "CC-BY-4.0", license_url = "https://creativecommons.org/licenses/by/4.0/",
    prepared_input_date = "2026-10-22", prepared_input_md5 = source_md5,
    preprocessing = c(
      "In-tissue spots matched to filtered Space Ranger counts; duplicate symbols collapsed before normalization.",
      "Retained genes detected in at least three spots and with nonzero normalized variance.",
      "log1p(10000 * counts / original spot library size), before gene filtering.",
      "Technical covariates and reference-informed composition estimates standardized; full-rank design retained.",
      "No additional spot QC, gene selection or covariate fitting in the public export."),
    export_changes = "Retained all expression genes and spots; included three complete Hallmark sets; removed duplicated arrays and environment-specific metadata.",
    gene_set_source = "MSigDB human Hallmark collection",
    gene_set_snapshot = "Frozen with the prepared input on 2026-10-22; exact MSigDB release was not recorded.",
    gene_set_license = "CC-BY-4.0",
    gene_set_license_url = "https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp",
    gene_set_copyright = "Copyright 2004-2025 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California"
  )
)
example <- passage_validate_bundle(example)
stopifnot(identical(example$Y, b$Y), identical(example$coords, b$coords),
          identical(example$X, b$X),
          identical(rownames(example$raw_counts), colnames(example$Y)),
          identical(colnames(example$raw_counts), rownames(example$Y)))
dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
saveRDS(example, args[2L], compress = "xz", version = 3)
check <- passage_validate_bundle(args[2L])
stopifnot(identical(check, example))
cat("Saved", nrow(check$Y), "spots and", ncol(check$Y), "genes\n")
print(tools::md5sum(args[2L]))
