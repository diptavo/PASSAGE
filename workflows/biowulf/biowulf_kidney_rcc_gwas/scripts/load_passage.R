passage_code <- Sys.getenv("PASSAGE_CODE", "")
source_roots <- unique(c(
  if (nzchar(passage_code)) passage_code else character(),
  file.path(root, "PASSAGE"),
  file.path(root, "SpaPath")
))
source_roots <- source_roots[dir.exists(file.path(source_roots, "R"))]

if (length(source_roots)) {
  r_dir <- file.path(source_roots[[1L]], "R")
  for (f in sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
  message("Loaded PASSAGE source from ", source_roots[[1L]])
} else if (requireNamespace("PASSAGE", quietly = TRUE)) {
  suppressPackageStartupMessages(library(PASSAGE))
  message("Loaded installed PASSAGE package")
} else {
  stop(
    "PASSAGE is unavailable. Set PASSAGE_CODE to the GitHub checkout ",
    "or install the PASSAGE R package."
  )
}
