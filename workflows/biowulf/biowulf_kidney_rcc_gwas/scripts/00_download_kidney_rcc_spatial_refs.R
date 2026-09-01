#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
setwd(root)
options(timeout = max(7200, getOption("timeout")))

dir.create(file.path(root, "data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "data", "spatial"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "refs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "refs", "gwas"), recursive = TRUE, showWarnings = FALSE)

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

zip_url <- "https://zenodo.org/records/14620362/files/TLS_VISIUM_USZ.zip?download=1"
zip_file <- file.path(root, "data", "raw", "TLS_VISIUM_USZ.zip")
message("Downloading Zenodo TLS kidney/lung cancer Visium archive")
z <- download_one(zip_url, zip_file)
if (!z$ok) stop("Failed to download spatial archive: ", z$message)

unzip_dir <- file.path(root, "data", "spatial", "TLS_VISIUM_USZ")
done <- file.path(unzip_dir, paste0(".", basename(zip_file), ".unzip_complete"))
if (!file.exists(done)) {
  dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
  message("Unzipping ", zip_file, " to ", unzip_dir)
  utils::unzip(zip_file, exdir = unzip_dir)
  writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), done)
} else {
  message("Spatial archive already unzipped: ", unzip_dir)
}

find_sample_dirs <- function(base_dir) {
  dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[basename(dirs) %in% c("KC1", "KC2", "KC3")]
  dirs <- dirs[file.exists(file.path(dirs, "filtered_feature_bc_matrix")) |
                 file.exists(file.path(dirs, "raw_feature_bc_matrix"))]
  dirs
}
sample_dirs <- sort(unique(find_sample_dirs(unzip_dir)))
if (!length(sample_dirs)) {
  stop("Could not locate KC1/KC2/KC3 Space Ranger output directories under ", unzip_dir)
}

spatial_manifest <- do.call(rbind, lapply(sample_dirs, function(d) {
  sample <- basename(d)
  data.frame(
    cancer = "kidney",
    disease = "renal_cell_carcinoma",
    spatial_sample = sample,
    source = "Zenodo 10x Visium kidney/lung cancer TLS dataset, record 14620362",
    source_url = zip_url,
    sample_dir = d,
    matrix_dir = if (dir.exists(file.path(d, "filtered_feature_bc_matrix"))) {
      file.path(d, "filtered_feature_bc_matrix")
    } else {
      file.path(d, "raw_feature_bc_matrix")
    },
    spatial_dir = file.path(d, "spatial"),
    annotation_file = {
      hits <- list.files(d, pattern = "TLS.*annotation|annotation.*[.]csv$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)
      if (length(hits)) hits[[1L]] else ""
    },
    stringsAsFactors = FALSE
  )
}))
write.csv(spatial_manifest, file.path(root, "data", "spatial_manifest.csv"), row.names = FALSE)

add_ref <- function(reference, cell_type, genes, source_note) {
  data.frame(cancer = "kidney", reference = reference, cell_type = cell_type,
             genes = paste(unique(toupper(genes)), collapse = ";"),
             source_note = source_note, stringsAsFactors = FALSE)
}

tme_tls <- list(
  malignant_ccrcc_hypoxia = c("CA9", "NDUFA4L2", "VEGFA", "EGLN3", "ANGPTL4", "SLC2A1", "ENO1"),
  malignant_epithelial = c("PAX8", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "AMACR"),
  t_nk = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1", "GZMB"),
  b_tls_plasma = c("MS4A1", "CD79A", "CD74", "CXCL13", "MZB1", "JCHAIN", "IGHG1"),
  myeloid = c("LST1", "TYROBP", "AIF1", "FCER1G", "C1QA", "C1QB", "SPP1"),
  fibroblast_stroma = c("COL1A1", "COL1A2", "DCN", "LUM", "TAGLN", "ACTA2"),
  endothelial = c("PECAM1", "VWF", "EMCN", "ENG", "KDR", "FLT1")
)
normal_kidney <- list(
  proximal_tubule = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "ALDOB", "AQP1"),
  loop_of_henle = c("UMOD", "SLC12A1", "CLDN16", "KCNJ1"),
  distal_tubule = c("SLC12A3", "CALB1", "TRPM6", "PVALB"),
  collecting_duct = c("AQP2", "AQP3", "KRT8", "KRT18", "ATP6V1B1", "SLC4A1"),
  podocyte = c("NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"),
  mesangial_pericyte = c("PDGFRB", "RGS5", "MCAM", "NOTCH3"),
  endothelial = c("PECAM1", "VWF", "EMCN", "ENG", "KDR"),
  immune = c("PTPRC", "CD3D", "MS4A1", "LST1", "NKG7"),
  stroma = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2")
)

rows <- list()
for (ct in names(tme_tls)) {
  rows[[length(rows) + 1L]] <- add_ref(
    "tumor_tme_tls_broad", ct, tme_tls[[ct]],
    "Broad RCC tumor/TME/TLS marker covariate panel for H3 adjustment."
  )
}
for (ct in names(normal_kidney)) {
  rows[[length(rows) + 1L]] <- add_ref(
    "normal_kidney_broad", ct, normal_kidney[[ct]],
    "Broad normal-kidney epithelial/stromal/immune marker covariate panel for H3 adjustment."
  )
}
marker_refs <- do.call(rbind, rows)
write.csv(marker_refs, file.path(root, "refs", "celltype_marker_references.csv"), row.names = FALSE)

gwas_glob <- Sys.getenv(
  "PASSAGE_GWAS_GLOB",
  "/data/Renal_GWAS_2022_exp/sumstats/meta_20230322/gwas_ssf/meta/meta.multianc.*.sumstats.tsv.bgz"
)
gwas_files <- Sys.glob(gwas_glob)
gwas_manifest <- data.frame(
  phenotype = sub("^meta[.]multianc[.]([^.]+)[.]sumstats[.]tsv[.]bgz$", "\\1", basename(gwas_files)),
  file = gwas_files,
  role = ifelse(grepl("[.]RCC[.]", basename(gwas_files)), "overall_RCC",
                ifelse(grepl("[.]CC2[.]", basename(gwas_files)), "clear_cell_RCC",
                       ifelse(grepl("[.]PRCC[.]", basename(gwas_files)), "papillary_RCC", "unknown"))),
  bytes = if (length(gwas_files)) file.info(gwas_files)$size else numeric(),
  stringsAsFactors = FALSE
)
gwas_manifest <- gwas_manifest[order(gwas_manifest$phenotype), , drop = FALSE]
write.csv(gwas_manifest, file.path(root, "refs", "gwas", "rcc_gwas_sumstats_manifest.csv"), row.names = FALSE)

writeLines(c(
  "# Kidney RCC PASSAGE Input Sources",
  "",
  paste0("- spatial_archive: ", zip_file),
  paste0("- spatial_samples: ", paste(spatial_manifest$spatial_sample, collapse = ", ")),
  paste0("- marker_references: ", paste(unique(marker_refs$reference), collapse = ", ")),
  paste0("- gwas_sumstats: ", paste(gwas_manifest$phenotype, collapse = ", "))
), file.path(root, "refs", "source_summary.md"))

message("Wrote spatial, cell-type marker, and GWAS manifests under ", root)
