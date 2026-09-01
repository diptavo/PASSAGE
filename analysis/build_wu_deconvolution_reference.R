# Build coarse cell-type pseudobulk signatures from Wu et al. GSE176078.
#
# The input matrix is Matrix Market genes x cells. This script streams the
# coordinate file and aggregates counts by cell type without materializing the
# full single-cell matrix in memory.
#
# Usage:
#   Rscript scripts/build_wu_deconvolution_reference.R

parse_args <- function(args) {
  cfg <- list(
    ref_dir = "data/reference/GSE176078/Wu_etal_2021_BRCA_scRNASeq",
    out_dir = "data/reference/GSE176078/processed",
    chunk_lines = 1000000L
  )
  for (arg in args) {
    if (grepl("^--ref-dir=", arg)) cfg$ref_dir <- sub("^--ref-dir=", "", arg)
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--chunk-lines=", arg)) cfg$chunk_lines <- as.integer(sub("^--chunk-lines=", "", arg))
  }
  cfg
}

coarse_cell_type <- function(x) {
  y <- rep(NA_character_, length(x))
  y[x == "Cancer Epithelial"] <- "malignant_epithelial"
  y[x == "Normal Epithelial"] <- "normal_epithelial"
  y[x == "CAFs"] <- "caf"
  y[x == "PVL"] <- "perivascular"
  y[x == "Endothelial"] <- "endothelial"
  y[x == "Myeloid"] <- "myeloid"
  y[x == "T-cells"] <- "t_nk"
  y[x == "B-cells"] <- "b_cell"
  y[x == "Plasmablasts"] <- "plasma_cell"
  y
}

read_mtx_header <- function(con) {
  first <- readLines(con, n = 1L)
  if (!grepl("^%%MatrixMarket", first)) {
    stop("Not a Matrix Market file")
  }
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) stop("Unexpected end of Matrix Market header")
    if (!startsWith(line, "%")) break
  }
  dims <- as.integer(strsplit(line, "\\s+")[[1]])
  if (length(dims) != 3L) stop("Could not parse Matrix Market dimensions")
  dims
}

aggregate_mtx_by_group <- function(mtx_file, n_genes, cell_group_index, n_groups, chunk_lines) {
  con <- file(mtx_file, open = "r")
  on.exit(close(con), add = TRUE)
  dims <- read_mtx_header(con)
  if (dims[1] != n_genes || dims[2] != length(cell_group_index)) {
    stop("Matrix dimensions do not match genes/barcodes: ",
         paste(dims, collapse = " x "))
  }
  message("Matrix dimensions: ", dims[1], " genes x ", dims[2], " cells; nnz=", dims[3])

  counts <- matrix(0, nrow = n_genes, ncol = n_groups)
  total_seen <- 0L
  repeat {
    x <- scan(
      con,
      what = list(i = integer(), j = integer(), v = numeric()),
      nlines = chunk_lines,
      quiet = TRUE
    )
    if (length(x$i) == 0L) break
    keep <- !is.na(cell_group_index[x$j])
    if (any(keep)) {
      idx <- x$i[keep] + (cell_group_index[x$j[keep]] - 1L) * n_genes
      sums <- tapply(x$v[keep], idx, sum)
      counts[as.integer(names(sums))] <- counts[as.integer(names(sums))] + as.numeric(sums)
    }
    total_seen <- total_seen + length(x$i)
    if (total_seen %% (10L * chunk_lines) < chunk_lines) {
      message("  streamed ", format(total_seen, big.mark = ","), " / ",
              format(dims[3], big.mark = ","), " nonzero entries")
    }
  }
  counts
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(cfg$ref_dir, "metadata.csv")
genes_file <- file.path(cfg$ref_dir, "count_matrix_genes.tsv")
barcodes_file <- file.path(cfg$ref_dir, "count_matrix_barcodes.tsv")
mtx_file <- file.path(cfg$ref_dir, "count_matrix_sparse.mtx")

message("Reading Wu metadata")
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
genes <- toupper(readLines(genes_file))
barcodes <- readLines(barcodes_file)
if (!identical(rownames(metadata), barcodes)) {
  stop("Metadata row names do not match matrix barcodes")
}

metadata$coarse_cell_type <- coarse_cell_type(metadata$celltype_major)
metadata <- metadata[!is.na(metadata$coarse_cell_type), , drop = FALSE]
cell_types <- sort(unique(metadata$coarse_cell_type))
message("Coarse cell types: ", paste(cell_types, collapse = ", "))
print(sort(table(metadata$coarse_cell_type), decreasing = TRUE))

cell_group <- rep(NA_integer_, length(barcodes))
names(cell_group) <- barcodes
cell_group[rownames(metadata)] <- match(metadata$coarse_cell_type, cell_types)

message("Streaming Matrix Market counts by coarse cell type")
counts_by_type <- aggregate_mtx_by_group(
  mtx_file = mtx_file,
  n_genes = length(genes),
  cell_group_index = cell_group,
  n_groups = length(cell_types),
  chunk_lines = cfg$chunk_lines
)
rownames(counts_by_type) <- genes
colnames(counts_by_type) <- cell_types

message("Collapsing duplicate gene symbols")
gene_groups <- factor(rownames(counts_by_type), levels = sort(unique(rownames(counts_by_type))))
collapsed <- rowsum(counts_by_type, group = gene_groups, reorder = FALSE)
cell_counts <- as.integer(table(factor(metadata$coarse_cell_type, levels = cell_types)))
names(cell_counts) <- cell_types
lib_sizes <- colSums(collapsed)

avg_counts_per_cell <- sweep(collapsed, 2L, pmax(cell_counts, 1L), "/")
cpm <- sweep(avg_counts_per_cell, 2L, pmax(colSums(avg_counts_per_cell), 1), "/") * 1e4
log_cpm <- log1p(cpm)

ref <- list(
  counts_by_type = collapsed,
  avg_counts_per_cell = avg_counts_per_cell,
  log_cpm = log_cpm,
  cell_counts = cell_counts,
  lib_sizes = lib_sizes,
  cell_types = cell_types,
  metadata_summary = data.frame(
    coarse_cell_type = cell_types,
    n_cells = as.integer(cell_counts),
    total_counts = as.numeric(lib_sizes),
    stringsAsFactors = FALSE
  ),
  source = "Wu_etal_2021_BRCA_scRNASeq_GSE176078"
)

saveRDS(ref, file.path(cfg$out_dir, "wu_brca_coarse_reference.rds"))
write.csv(ref$metadata_summary, file.path(cfg$out_dir, "wu_brca_coarse_reference_summary.csv"), row.names = FALSE)
write.csv(ref$log_cpm, file.path(cfg$out_dir, "wu_brca_coarse_logcpm.csv"))
message("Wrote reference to ", cfg$out_dir)
