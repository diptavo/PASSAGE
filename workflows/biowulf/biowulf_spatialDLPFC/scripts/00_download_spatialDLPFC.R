#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
data_dir <- file.path(root, "data", "spatialDLPFC")
ref_dir <- file.path(root, "refs", "spatialDLPFC_snRNAseq")
cache_dir <- file.path(root, "cache", "spatialLIBD")
hub_cache_dir <- file.path(root, "cache", "ExperimentHub")
r_user_cache_dir <- file.path(root, "cache", "R_user_cache")
bfc_cache_dir <- file.path(root, "cache", "BiocFileCache")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hub_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_user_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bfc_cache_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  EXPERIMENT_HUB_CACHE = hub_cache_dir,
  ANNOTATION_HUB_CACHE = hub_cache_dir,
  R_USER_CACHE_DIR = r_user_cache_dir
)
options(
  ExperimentHub.Cache = hub_cache_dir,
  AnnotationHub.Cache = hub_cache_dir,
  BiocFileCache.cache = bfc_cache_dir
)

message("Root: ", root)
message("Data dir: ", data_dir)
message("Reference dir: ", ref_dir)
message("ExperimentHub cache dir: ", hub_cache_dir)
message("R user cache dir: ", r_user_cache_dir)
message("BiocFileCache dir: ", bfc_cache_dir)

ensure_pkg <- function(pkg, bioc = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (bioc) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

for (pkg in c("spatialLIBD", "SpatialExperiment", "SingleCellExperiment",
              "SummarizedExperiment", "HDF5Array", "BiocFileCache",
              "ExperimentHub")) {
  ensure_pkg(pkg, bioc = TRUE)
}

suppressPackageStartupMessages({
  library(spatialLIBD)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(S4Vectors)
})

drop_spatial_images <- function(spe) {
  out <- tryCatch({
    SpatialExperiment::imgData(spe) <- S4Vectors::DataFrame(
      sample_id = character(),
      image_id = character(),
      data = I(list()),
      scaleFactor = numeric()
    )
    spe
  }, error = function(e) {
    message("Could not reset imgData with empty DataFrame: ", conditionMessage(e))
    tryCatch({
      SpatialExperiment::imgData(spe) <- NULL
      spe
    }, error = function(e2) {
      message("Could not drop imgData: ", conditionMessage(e2))
      spe
    })
  })
  out
}

spe_rds <- file.path(data_dir, "spatialDLPFC_Visium_spe.rds")
spe <- NULL
if (file.exists(spe_rds)) {
  spe <- tryCatch(readRDS(spe_rds), error = function(e) {
    message("Existing Visium RDS is not readable, replacing it: ", conditionMessage(e))
    unlink(spe_rds)
    NULL
  })
}
if (is.null(spe)) {
  message("Downloading spatialDLPFC Visium SpatialExperiment")
  spe <- spatialLIBD::fetch_data(
    type = "spatialDLPFC_Visium",
    destdir = cache_dir
  )
  spe <- drop_spatial_images(spe)
  saveRDS(spe, spe_rds, compress = FALSE)
} else {
  message("Visium RDS already exists: ", spe_rds)
}

visium_summary <- data.frame(
  object = "spatialDLPFC_Visium",
  genes = nrow(spe),
  spots = ncol(spe),
  assays = paste0(assayNames(spe), collapse = ","),
  coldata_columns = ncol(colData(spe)),
  stringsAsFactors = FALSE
)
write.csv(visium_summary, file.path(data_dir, "spatialDLPFC_Visium_summary.csv"),
          row.names = FALSE)

sample_col <- intersect(c("sample_id", "sample_name", "sample", "array"),
                        colnames(colData(spe)))[1L]
if (!is.na(sample_col)) {
  sample_tab <- as.data.frame(table(colData(spe)[[sample_col]]))
  colnames(sample_tab) <- c(sample_col, "n_spots")
  write.csv(sample_tab, file.path(data_dir, "spatialDLPFC_Visium_samples.csv"),
            row.names = FALSE)
}

sn_zip_record <- file.path(ref_dir, "spatialDLPFC_snRNAseq_fetch_path.txt")
if (!file.exists(sn_zip_record)) {
  message("Downloading spatialDLPFC snRNA-seq reference")
  sn_zip <- spatialLIBD::fetch_data(
    type = "spatialDLPFC_snRNAseq",
    destdir = cache_dir
  )
  writeLines(sn_zip, sn_zip_record)
  unzip(sn_zip, exdir = ref_dir)
} else {
  message("snRNA-seq reference fetch record already exists: ", sn_zip_record)
}

sce_dir <- file.path(ref_dir, "sce_DLPFC_annotated")
if (dir.exists(sce_dir)) {
  message("Loading HDF5-backed snRNA-seq reference")
  sce <- HDF5Array::loadHDF5SummarizedExperiment(sce_dir)
  sn_summary <- data.frame(
    object = "spatialDLPFC_snRNAseq",
    genes = nrow(sce),
    cells = ncol(sce),
    assays = paste0(assayNames(sce), collapse = ","),
    coldata_columns = ncol(colData(sce)),
    stringsAsFactors = FALSE
  )
  write.csv(sn_summary, file.path(ref_dir, "spatialDLPFC_snRNAseq_summary.csv"),
            row.names = FALSE)
  cd <- as.data.frame(colData(sce))
  cell_cols <- grep("cell|type|cluster|annot|layer", colnames(cd),
                    ignore.case = TRUE, value = TRUE)
  write.csv(cd[, cell_cols, drop = FALSE],
            file.path(ref_dir, "spatialDLPFC_snRNAseq_cell_annotations.csv"),
            row.names = TRUE)
}

message("Download completed")
