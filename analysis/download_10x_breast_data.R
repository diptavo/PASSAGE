# Download the two public 10x Visium breast cancer datasets used by the
# SpaPath Hallmark analysis.
#
# Usage:
#   cd /path/to/PASSAGE
#   Rscript scripts/download_10x_breast_data.R

files <- data.frame(
  dataset = c(
    "Visium_FFPE_Human_Breast_Cancer",
    "Visium_FFPE_Human_Breast_Cancer",
    "Visium_FFPE_Human_Breast_Cancer",
    "V1_Breast_Cancer_Block_A_Section_1",
    "V1_Breast_Cancer_Block_A_Section_1"
  ),
  file = c(
    "filtered_feature_bc_matrix.h5",
    "spatial.tar.gz",
    "pathologist_annotations.png",
    "filtered_feature_bc_matrix.h5",
    "spatial.tar.gz"
  ),
  url = c(
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(files))) {
  out_dir <- file.path("data", "raw", files$dataset[i])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, files$file[i])
  if (!file.exists(out_file)) {
    message("Downloading ", files$url[i])
    download.file(files$url[i], out_file, mode = "wb", quiet = FALSE)
  } else {
    message("Already exists: ", out_file)
  }

  if (grepl("spatial[.]tar[.]gz$", out_file)) {
    spatial_dir <- file.path(out_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
      message("Extracting ", out_file)
      utils::untar(out_file, exdir = out_dir)
    }
  }
}
