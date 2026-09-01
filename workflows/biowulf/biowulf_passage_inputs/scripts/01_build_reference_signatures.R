#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "/data/DCEG_Dutta/PASSAGE_production_inputs_20260825"
source(file.path(root, "scripts", "reference_common.R"))

dir.create(file.path(root, "refs", "derived"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "metadata"), recursive = TRUE, showWarnings = FALSE)

build_breast <- function() {
  ref_dir <- file.path(root, "refs", "breast", "GSE176078", "Wu_etal_2021_BRCA_scRNASeq")
  counts <- readMM(file.path(ref_dir, "count_matrix_sparse.mtx"))
  genes <- fread(file.path(ref_dir, "count_matrix_genes.tsv"), header = FALSE, data.table = FALSE)[[1L]]
  barcodes <- fread(file.path(ref_dir, "count_matrix_barcodes.tsv"), header = FALSE, data.table = FALSE)[[1L]]
  if (nrow(counts) != length(genes) && ncol(counts) == length(genes)) counts <- t(counts)
  stopifnot(nrow(counts) == length(genes), ncol(counts) == length(barcodes))
  rownames(counts) <- clean_symbol(genes)
  colnames(counts) <- barcodes
  counts <- collapse_sparse_rows(counts, rownames(counts))

  meta <- read.csv(file.path(ref_dir, "metadata.csv"), row.names = 1L, stringsAsFactors = FALSE, check.names = FALSE)
  idx <- match(colnames(counts), rownames(meta))
  stopifnot(mean(!is.na(idx)) > 0.95)
  counts <- counts[, !is.na(idx), drop = FALSE]
  meta <- meta[idx[!is.na(idx)], , drop = FALSE]
  major_col <- intersect(c("celltype_major", "celltype_subset", "celltype_minor"), names(meta))[[1L]]
  cell_type <- map_breast_type(meta[[major_col]])
  result <- group_pseudobulk_signature(counts, as.character(meta$orig.ident), cell_type)
  audit <- data.frame(
    original_label = as.character(meta[[major_col]]),
    broad_cell_type = cell_type,
    stringsAsFactors = FALSE
  )
  out <- list(
    cohort = "breast",
    reference_id = "GSE176078",
    source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
    method = "donor-balanced pseudobulk CPM reference; data-driven and curated marker genes",
    signature = result$signature,
    group_table = result$group_table,
    cell_counts = result$cell_counts,
    label_mapping = as.data.frame(table(audit), stringsAsFactors = FALSE)
  )
  saveRDS(out, file.path(root, "refs", "derived", "breast_GSE176078_signature.rds"), compress = "gzip")
  out
}

build_kidney <- function() {
  ref_dir <- file.path(root, "refs", "kidney", "GSE224630")
  data_file <- file.path(ref_dir, "GSE224630_normalized_integrated.data.tsv.gz")
  meta <- read_table_any(file.path(ref_dir, "GSE224630_overall_metadata.tsv.gz"))
  names(meta) <- make.names(names(meta))
  stopifnot(all(c("barcode", "patient", "state", "all.cell.association") %in% names(meta)))
  cell_type <- map_kidney_type(meta$all.cell.association, meta$state)

  marker_file <- file.path(ref_dir, "kidney_signature_marker_genes.txt")
  marker_genes <- sort(unique(unlist(canonical_markers$kidney, use.names = FALSE)))
  writeLines(marker_genes, marker_file)
  header <- strsplit(readLines(gzfile(data_file), n = 1L), "\t", fixed = TRUE)[[1L]]
  awk_program <- "NR==FNR {keep[toupper($1)]=1; next} (toupper($1) in keep)"
  cmd <- paste("gzip -dc", shQuote(data_file), "| awk -F '\\t'", shQuote(awk_program), shQuote(marker_file), "-")
  selected <- fread(cmd = cmd, header = FALSE, data.table = FALSE, showProgress = FALSE)
  stopifnot(ncol(selected) == length(header) + 1L)
  genes <- clean_symbol(selected[[1L]])
  values <- as.matrix(selected[, -1L, drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- genes
  colnames(values) <- header

  idx <- match(colnames(values), meta$barcode)
  stopifnot(mean(!is.na(idx)) > 0.90)
  values <- values[, !is.na(idx), drop = FALSE]
  meta2 <- meta[idx[!is.na(idx)], , drop = FALSE]
  cell_type2 <- map_kidney_type(meta2$all.cell.association, meta2$state)
  result <- donor_balanced_dense_signature(values, as.character(meta2$patient), cell_type2)
  audit <- data.frame(
    original_label = as.character(meta2$all.cell.association),
    state = as.character(meta2$state),
    broad_cell_type = cell_type2,
    stringsAsFactors = FALSE
  )
  out <- list(
    cohort = "kidney",
    reference_id = "GSE224630",
    source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630",
    method = "donor-balanced annotated-cell reference restricted to curated kidney and ccRCC marker genes",
    signature = result$signature,
    group_table = result$group_table,
    cell_counts = result$cell_counts,
    label_mapping = as.data.frame(table(audit), stringsAsFactors = FALSE)
  )
  saveRDS(out, file.path(root, "refs", "derived", "kidney_GSE224630_signature.rds"), compress = "gzip")
  out
}

breast <- build_breast()
kidney <- build_kidney()
summary <- do.call(rbind, lapply(list(breast, kidney), function(x) data.frame(
  cohort = x$cohort,
  reference_id = x$reference_id,
  n_genes = nrow(x$signature),
  n_cell_types = ncol(x$signature),
  cell_types = paste(colnames(x$signature), collapse = ";"),
  n_donors = length(unique(x$group_table$donor)),
  stringsAsFactors = FALSE
)))
write.csv(summary, file.path(root, "metadata", "reference_signature_summary.csv"), row.names = FALSE)
date <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
writeLines(date, file.path(root, "metadata", "reference_signatures.complete"))
