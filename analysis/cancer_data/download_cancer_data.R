# Download selected complete samples, or all samples, with integrity checks.
download_cancer_data <- function(samples = NULL, dest = "cancer_data", ref = "main") {
  base <- paste0("https://raw.githubusercontent.com/diptavo/PASSAGE/", ref,
                 "/analysis/cancer_data/")
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  manifest_file <- tempfile(fileext = ".csv")
  on.exit(unlink(manifest_file), add = TRUE)
  download.file(paste0(base, "manifest.csv"), manifest_file, mode = "wb", quiet = TRUE)
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
  required <- c("sample", "relative_path", "md5", "n_spots", "n_genes")
  if (!all(required %in% names(manifest)) || anyDuplicated(manifest$sample)) {
    stop("Invalid dataset manifest")
  }
  if (!is.null(samples)) {
    missing <- setdiff(samples, manifest$sample)
    if (length(missing)) stop("Unknown sample(s): ", paste(missing, collapse = ", "))
    manifest <- manifest[match(unique(samples), manifest$sample), , drop = FALSE]
  }
  if (!nrow(manifest)) stop("No samples selected")
  if (any(!grepl("^[a-z]+/[A-Za-z0-9_]+[.]rds$", manifest$relative_path))) {
    stop("Unexpected dataset path in manifest")
  }
  manifest$local_file <- file.path(dest, manifest$relative_path)
  for (i in seq_len(nrow(manifest))) {
    path <- manifest$local_file[i]
    if (file.exists(path)) {
      if (unname(tools::md5sum(path)) != manifest$md5[i]) {
        stop("Existing file has a different checksum; choose a new destination: ", path)
      }
      message("Already verified: ", manifest$sample[i])
      next
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    part <- tempfile(pattern = "download_", tmpdir = dirname(path))
    on.exit(unlink(part), add = TRUE)
    message("Downloading: ", manifest$sample[i])
    download.file(paste0(base, manifest$relative_path[i]), part, mode = "wb")
    if (unname(tools::md5sum(part)) != manifest$md5[i]) stop("Checksum mismatch: ", path)
    if (!file.rename(part, path)) stop("Could not save downloaded file: ", path)
  }
  invisible(manifest)
}
