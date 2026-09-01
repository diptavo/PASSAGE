#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)
options(timeout = max(3600, getOption("timeout")))

dir.create(file.path(root, "data", "spatial"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "refs", "raw_scRNA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "refs"), recursive = TRUE, showWarnings = FALSE)

spatial <- data.frame(
  cancer = c("breast", "cervical", "prostate", "lung"),
  spatial_sample = c("Visium_FFPE_Human_Breast_Cancer",
                     "Visium_FFPE_Human_Cervical_Cancer",
                     "Visium_FFPE_Human_Prostate_IF",
                     "GSM8855706_NSCLC_P4"),
  source = c("10x Genomics public Visium FFPE breast cancer",
             "10x Genomics public Visium FFPE cervical cancer",
             "10x Genomics public Visium FFPE prostate cancer",
             "GEO GSE292299/GSM8855706 NSCLC patient 4 Visium"),
  matrix_url = c(
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8855nnn/GSM8855706/suppl/GSM8855706_NSCLC_P4_filtered_feature_bc_matrix.tar.gz"
  ),
  spatial_url = c(
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_spatial.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8855nnn/GSM8855706/suppl/GSM8855706_NSCLC_P4_spatial.tar.gz"
  ),
  stringsAsFactors = FALSE
)

raw_refs <- data.frame(
  cancer = c("breast", "cervical", "prostate", "lung"),
  reference_id = c("GSE176078_Wu_BRCA_scRNA",
                   "GSE168652_cervical_scRNA_RAW",
                   "GSE176031_prostate_scRNA_RAW",
                   "GSE131907_lung_cancer_scRNA_RAW"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168652/suppl/GSE168652_RAW.tar",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176031/suppl/GSE176031_RAW.tar",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_RAW.tar"
  ),
  role = c("tumor immune/stromal/malignant breast cancer reference",
           "tumor immune/stromal/malignant cervical cancer reference",
           "tumor immune/stromal/malignant prostate cancer reference",
           "tumor immune/stromal/malignant lung cancer reference"),
  stringsAsFactors = FALSE
)

download_one <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(list(ok = TRUE, message = "exists", bytes = file.info(dest)$size))
  }
  msg <- tryCatch({
    utils::download.file(url, dest, mode = "wb", method = "libcurl", quiet = FALSE)
    "downloaded"
  }, error = function(e) conditionMessage(e))
  ok <- file.exists(dest) && file.info(dest)$size > 0
  list(ok = ok, message = msg, bytes = if (ok) file.info(dest)$size else 0)
}

untar_if_needed <- function(tarfile, exdir) {
  done <- file.path(exdir, paste0(".", basename(tarfile), ".untar_complete"))
  if (file.exists(done)) return(TRUE)
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    utils::untar(tarfile, exdir = exdir)
    writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), done)
    TRUE
  }, error = function(e) {
    message("untar failed for ", tarfile, ": ", conditionMessage(e))
    FALSE
  })
  ok
}

spatial_status <- vector("list", nrow(spatial))
for (ii in seq_len(nrow(spatial))) {
  row <- spatial[ii, ]
  sample_dir <- file.path(root, "data", "spatial", row$cancer, row$spatial_sample)
  dl_dir <- file.path(sample_dir, "downloads")
  outs_dir <- file.path(sample_dir, "outs")
  matrix_tar <- file.path(dl_dir, basename(row$matrix_url))
  spatial_tar <- file.path(dl_dir, basename(row$spatial_url))
  message("Downloading spatial data for ", row$cancer, " / ", row$spatial_sample)
  m <- download_one(row$matrix_url, matrix_tar)
  s <- download_one(row$spatial_url, spatial_tar)
  matrix_untar <- if (m$ok) untar_if_needed(matrix_tar, outs_dir) else FALSE
  spatial_untar <- if (s$ok) untar_if_needed(spatial_tar, outs_dir) else FALSE
  spatial_status[[ii]] <- cbind(
    row,
    matrix_file = matrix_tar,
    spatial_file = spatial_tar,
    sample_dir = sample_dir,
    outs_dir = outs_dir,
    matrix_ok = m$ok,
    spatial_ok = s$ok,
    matrix_untar_ok = matrix_untar,
    spatial_untar_ok = spatial_untar,
    matrix_message = m$message,
    spatial_message = s$message,
    stringsAsFactors = FALSE
  )
}
spatial_manifest <- do.call(rbind, spatial_status)
write.csv(spatial_manifest, file.path(root, "data", "spatial_manifest.csv"), row.names = FALSE)

ref_status <- vector("list", nrow(raw_refs))
for (ii in seq_len(nrow(raw_refs))) {
  row <- raw_refs[ii, ]
  dest <- file.path(root, "refs", "raw_scRNA", row$cancer, basename(row$url))
  message("Staging raw scRNA reference ", row$reference_id)
  x <- download_one(row$url, dest)
  ref_status[[ii]] <- cbind(row, file = dest, ok = x$ok, bytes = x$bytes,
                            message = x$message, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, ref_status),
          file.path(root, "refs", "raw_scRNA", "download_status.csv"),
          row.names = FALSE)

add_ref <- function(cancer, reference, cell_type, genes, source_note) {
  data.frame(cancer = cancer, reference = reference, cell_type = cell_type,
             genes = paste(unique(toupper(genes)), collapse = ";"),
             source_note = source_note, stringsAsFactors = FALSE)
}

rows <- list()
common_tme <- list(
  t_nk = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"),
  b_plasma = c("MS4A1", "CD79A", "CD74", "MZB1", "JCHAIN", "IGHG1"),
  myeloid = c("LST1", "TYROBP", "AIF1", "FCER1G", "C1QA", "C1QB", "SPP1"),
  fibroblast_caf = c("COL1A1", "COL1A2", "DCN", "LUM", "TAGLN", "ACTA2"),
  endothelial = c("PECAM1", "VWF", "ENG", "KDR"),
  epithelial_malignant = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MKI67", "TOP2A")
)
normal_shared <- list(
  immune = c("PTPRC", "CD3D", "MS4A1", "LST1", "NKG7"),
  stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "PECAM1", "VWF"),
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7")
)
tissue <- list(
  breast = list(
    tumor_tme_broad = c(common_tme, list(
      luminal = c("ESR1", "PGR", "FOXA1", "GATA3", "KRT8", "KRT18"),
      basal_myoepithelial = c("KRT5", "KRT14", "KRT17", "TP63", "ACTA2", "MYH11"),
      adipocyte = c("ADIPOQ", "PLIN1", "LPL", "FABP4")
    )),
    normal_tissue_broad = c(normal_shared, list(
      luminal_epithelial = c("KRT8", "KRT18", "KRT19", "EPCAM", "ESR1"),
      basal_myoepithelial = c("KRT5", "KRT14", "TP63", "ACTA2"),
      adipocyte = c("ADIPOQ", "PLIN1", "LPL", "FABP4")
    ))
  ),
  cervical = list(
    tumor_tme_broad = c(common_tme, list(
      squamous_epithelial = c("KRT5", "KRT14", "KRT17", "TP63", "SFN"),
      glandular_epithelial = c("KRT7", "KRT8", "KRT18", "KRT19", "MUC1"),
      smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1")
    )),
    normal_tissue_broad = c(normal_shared, list(
      squamous_epithelial = c("KRT5", "KRT14", "KRT17", "TP63"),
      glandular_epithelial = c("KRT7", "KRT8", "KRT18", "KRT19", "MUC1"),
      smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1")
    ))
  ),
  prostate = list(
    tumor_tme_broad = c(common_tme, list(
      luminal_prostate = c("KLK3", "ACPP", "NKX3-1", "KRT8", "KRT18", "AR"),
      basal_prostate = c("KRT5", "KRT14", "TP63"),
      club_secretory = c("SCGB3A1", "PIGR", "MUC1")
    )),
    normal_tissue_broad = c(normal_shared, list(
      luminal_prostate = c("KLK3", "ACPP", "NKX3-1", "KRT8", "KRT18"),
      basal_prostate = c("KRT5", "KRT14", "TP63"),
      smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1")
    ))
  ),
  lung = list(
    tumor_tme_broad = c(common_tme, list(
      alveolar_type2 = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "NAPSA"),
      club_ciliated = c("SCGB1A1", "SCGB3A1", "FOXJ1", "PIFO"),
      basal_squamous = c("KRT5", "KRT14", "KRT17", "TP63")
    )),
    normal_tissue_broad = c(normal_shared, list(
      alveolar_type2 = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "NAPSA"),
      club_ciliated = c("SCGB1A1", "SCGB3A1", "FOXJ1", "PIFO"),
      basal = c("KRT5", "KRT14", "TP63")
    ))
  )
)
for (cc in names(tissue)) {
  for (rr in names(tissue[[cc]])) {
    marker_sets <- tissue[[cc]][[rr]]
    for (ct in names(marker_sets)) {
      rows[[length(rows) + 1L]] <- add_ref(
        cc, rr, ct, marker_sets[[ct]],
        "Broad marker-derived covariate reference; raw scRNA archives are staged separately for deconvolution follow-up."
      )
    }
  }
}
marker_refs <- do.call(rbind, rows)
write.csv(marker_refs, file.path(root, "refs", "celltype_marker_references.csv"),
          row.names = FALSE)

message("Wrote spatial manifest and marker references under ", root)
