# Real-data smoke comparison for PASSAGE anti-double-dipping modes on the two
# public 10x Visium breast cancer datasets.
#
# This is intentionally not the full production analysis. It writes compact CSV
# and Markdown summaries only, with small permutation counts suitable for local
# feasibility checks.
#
# Usage:
#   Rscript scripts/run_10x_breast_antidip_smoke.R \
#     --n-perm=19 --n-perm-refit=5 --max-pathways=8 \
#     --out-dir=results/passage_10x_breast_antidip_smoke

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
})

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

parse_args <- function(args) {
  cfg <- list(
    n_perm = 19L,
    n_perm_refit = 5L,
    K = 4L,
    m = 10L,
    n_folds = 4L,
    seed = 20260730L,
    max_pathways = 8L,
    max_refit_pathways = 3L,
    out_dir = file.path("results", "passage_10x_breast_antidip_smoke"),
    covariate_mode = "technical_spatial_quadratic",
    gmt = file.path("msigdb", "h.all.v2023.2.Hs.symbols.gmt"),
    datasets = c("Visium_FFPE_Human_Breast_Cancer", "V1_Breast_Cancer_Block_A_Section_1"),
    modes = c("none", "spot_crossfit", "pathway_holdout", "refit_null")
  )
  for (arg in args) {
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--n-perm-refit=", arg)) cfg$n_perm_refit <- as.integer(sub("^--n-perm-refit=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--max-refit-pathways=", arg)) cfg$max_refit_pathways <- as.integer(sub("^--max-refit-pathways=", "", arg))
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--covariate-mode=", arg)) cfg$covariate_mode <- sub("^--covariate-mode=", "", arg)
    if (grepl("^--gmt=", arg)) cfg$gmt <- sub("^--gmt=", "", arg)
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--modes=", arg)) cfg$modes <- strsplit(sub("^--modes=", "", arg), ",", fixed = TRUE)[[1]]
  }
  cfg
}

dataset_configs <- list(
  Visium_FFPE_Human_Breast_Cancer = list(
    label = "10x Visium FFPE Human Breast Cancer DCIS/Invasive",
    root = file.path("data", "raw", "Visium_FFPE_Human_Breast_Cancer")
  ),
  V1_Breast_Cancer_Block_A_Section_1 = list(
    label = "10x Visium Breast Cancer Block A Section 1",
    root = file.path("data", "raw", "V1_Breast_Cancer_Block_A_Section_1")
  )
)

read_gmt_pathways <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    alt <- file.path("..", "PASSAGE_run_bundle", "SpaPath", gmt_path)
    if (file.exists(alt)) gmt_path <- alt
  }
  if (!file.exists(gmt_path)) {
    stop("GMT file not found: ", gmt_path)
  }
  lines <- readLines(gmt_path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) next
    out[[fields[[1L]]]] <- toupper(unique(fields[-c(1L, 2L)]))
  }
  out
}

read_10x_h5 <- function(h5_path) {
  data <- rhdf5::h5read(h5_path, "/matrix/data")
  indices <- rhdf5::h5read(h5_path, "/matrix/indices")
  indptr <- rhdf5::h5read(h5_path, "/matrix/indptr")
  shape <- rhdf5::h5read(h5_path, "/matrix/shape")
  barcodes <- rhdf5::h5read(h5_path, "/matrix/barcodes")
  genes <- rhdf5::h5read(h5_path, "/matrix/features/name")
  mat <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = as.integer(shape)
  )
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  as(mat, "dgCMatrix")
}

read_tissue_positions <- function(spatial_dir) {
  f <- file.path(spatial_dir, "tissue_positions.csv")
  header <- TRUE
  if (!file.exists(f)) {
    f <- file.path(spatial_dir, "tissue_positions_list.csv")
    header <- FALSE
  }
  if (!file.exists(f)) stop("Could not find tissue positions file in ", spatial_dir)
  pos <- read.csv(f, header = header, stringsAsFactors = FALSE)
  if (!header) {
    colnames(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }
  pos
}

collapse_counts_by_symbol <- function(counts) {
  symbols <- toupper(rownames(counts))
  ok <- !is.na(symbols) & nzchar(symbols)
  counts <- counts[ok, , drop = FALSE]
  symbols <- symbols[ok]
  groups <- factor(symbols, levels = sort(unique(symbols)))
  G <- Matrix::sparse.model.matrix(~ 0 + groups)
  out <- t(G) %*% counts
  rownames(out) <- levels(groups)
  colnames(out) <- colnames(counts)
  as(out, "dgCMatrix")
}

prepare_dataset <- function(config, pathways_input, covariate_mode = "technical_spatial_quadratic",
                            min_detected_fraction = 0.01, scale_factor = 1e4) {
  counts_raw <- read_10x_h5(file.path(config$root, "filtered_feature_bc_matrix.h5"))
  pos <- read_tissue_positions(file.path(config$root, "spatial"))
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(counts_raw), pos$barcode)
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]
  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)
  counts <- collapse_counts_by_symbol(counts_raw)
  pathway_genes <- sort(unique(unlist(pathways_input, use.names = FALSE)))
  counts <- counts[intersect(pathway_genes, rownames(counts)), , drop = FALSE]
  min_spots <- max(20, ceiling(min_detected_fraction * ncol(counts)))
  counts <- counts[Matrix::rowSums(counts > 0) >= min_spots, , drop = FALSE]
  norm <- counts %*% Matrix::Diagonal(x = scale_factor / pmax(lib_size, 1))
  norm@x <- log1p(norm@x)
  Y <- t(as.matrix(norm))
  Y <- Y[, apply(Y, 2L, stats::var) > 1e-8, drop = FALSE]
  coords <- as.matrix(pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  coords <- scale(coords)
  coords <- sweep(coords, 2L, apply(coords, 2L, min), "-")
  coords <- sweep(coords, 2L, apply(coords, 2L, max), "/")
  X <- cbind(
    log_umi = as.numeric(scale(log1p(lib_size))),
    log_detected = as.numeric(scale(log1p(n_detected)))
  )
  if (covariate_mode %in% c("technical_spatial_linear", "technical_spatial_quadratic")) {
    sx <- as.numeric(scale(coords[, 1L]))
    sy <- as.numeric(scale(coords[, 2L]))
    X <- cbind(X, spatial_x = sx, spatial_y = sy)
  }
  if (covariate_mode == "technical_spatial_quadratic") {
    sx <- X[, "spatial_x"]
    sy <- X[, "spatial_y"]
    X <- cbind(
      X,
      spatial_x2 = as.numeric(scale(sx^2)),
      spatial_y2 = as.numeric(scale(sy^2)),
      spatial_xy = as.numeric(scale(sx * sy))
    )
  }
  rownames(X) <- common
  pathways <- lapply(pathways_input, intersect, y = colnames(Y))
  pathways <- pathways[lengths(pathways) >= 5]
  list(Y = Y, coords = coords, X = X, pathways = pathways, n_spots = nrow(Y), n_features = ncol(Y))
}

select_pathways <- function(pathways, max_pathways) {
  preferred <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_GLYCOLYSIS"
  )
  keep <- intersect(preferred, names(pathways))
  keep <- c(keep, setdiff(names(pathways), keep))
  pathways[keep[seq_len(min(length(keep), max_pathways))]]
}

run_mode <- function(dat, dataset_id, mode, cfg) {
  pw <- dat$pathways
  if (mode == "refit_null") {
    pw <- pw[seq_len(min(length(pw), cfg$max_refit_pathways))]
  }
  range_grid <- passage_default_range_grid(dat$coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
  fit_args <- list(
    K = cfg$K,
    m = cfg$m,
    range_grid = range_grid,
    kernel = "matern32",
    ordering = "none",
    verbose = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- if (mode == "none") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = pw,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = cfg$n_perm,
      seed = cfg$seed, verbose = FALSE
    )
  } else if (mode == "spot_crossfit") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = pw,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = cfg$n_perm,
      seed = cfg$seed, verbose = FALSE, anti_dip = "spot_crossfit",
      anti_dip_args = list(n_folds = cfg$n_folds, seed = cfg$seed, fit_args = fit_args)
    )
  } else if (mode == "pathway_holdout") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = pw,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = cfg$n_perm,
      seed = cfg$seed, verbose = FALSE, anti_dip = "pathway_holdout",
      anti_dip_args = list(fit_args = fit_args)
    )
  } else if (mode == "refit_null") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = pw,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "moment", n_perm = cfg$n_perm_refit,
      seed = cfg$seed, verbose = FALSE, anti_dip = "refit_null",
      anti_dip_args = list(fit_args = fit_args, n_perm = cfg$n_perm_refit)
    )
  } else {
    stop("Unknown mode: ", mode)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  tbl <- fit$summary
  if ("spasset_genes_H1" %in% names(tbl) && is.list(tbl$spasset_genes_H1)) {
    tbl$spasset_genes_H1 <- vapply(tbl$spasset_genes_H1, function(x) {
      paste(as.character(x), collapse = ";")
    }, character(1))
  }
  tbl$dataset <- dataset_id
  tbl$mode <- mode
  tbl$total_mode_elapsed_sec <- elapsed
  tbl$n_spots <- dat$n_spots
  tbl$n_features <- dat$n_features
  tbl$n_perm <- if (mode == "refit_null") cfg$n_perm_refit else cfg$n_perm
  tbl
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
pathways_all <- read_gmt_pathways(cfg$gmt)
rows <- list()
ii <- 0L
for (dataset_id in cfg$datasets) {
  message("\n=== Preparing ", dataset_id, " ===")
  dat <- prepare_dataset(dataset_configs[[dataset_id]], pathways_all, cfg$covariate_mode)
  dat$pathways <- select_pathways(dat$pathways, cfg$max_pathways)
  message("Prepared ", dat$n_spots, " spots, ", dat$n_features, " genes; selected ",
          length(dat$pathways), " pathways")
  for (mode in cfg$modes) {
    message("Running mode=", mode)
    ii <- ii + 1L
    rows[[ii]] <- run_mode(dat, dataset_id, mode, cfg)
  }
}
out <- do.call(rbind, rows)
out$fdr_by_dataset_mode <- ave(out$p_H1, out$dataset, out$mode, FUN = function(p) stats::p.adjust(p, "BH"))
out <- out[order(out$dataset, out$mode, out$p_H1), , drop = FALSE]
csv <- file.path(cfg$out_dir, "breast_antidip_pathway_results.csv")
write.csv(out, csv, row.names = FALSE)

summary_tbl <- aggregate(
  cbind(p_H1, R2_cca, PSVS_range, mean_propSV, total_mode_elapsed_sec) ~ dataset + mode,
  data = out,
  FUN = function(x) c(median = stats::median(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
)
write.csv(summary_tbl, file.path(cfg$out_dir, "breast_antidip_mode_summary.csv"), row.names = FALSE)

top_cols <- c("dataset", "mode", "pathway", "pathway_size", "p_H1", "fdr_by_dataset_mode",
              "p_H1_moment", "R2_cca", "PSVS_range", "mean_propSV", "total_mode_elapsed_sec")
top <- out[, intersect(top_cols, names(out)), drop = FALSE]
md <- c(
  "# PASSAGE Breast Cancer Anti-Dip Smoke",
  "",
  paste0("- Datasets: ", paste(cfg$datasets, collapse = ", ")),
  paste0("- Modes: ", paste(cfg$modes, collapse = ", ")),
  paste0("- Pathways per non-refit mode: ", cfg$max_pathways),
  paste0("- Refit-null pathways per dataset: ", cfg$max_refit_pathways),
  paste0("- Permutations: ", cfg$n_perm, " non-refit; ", cfg$n_perm_refit, " refit-null"),
  paste0("- K: ", cfg$K, "; m: ", cfg$m, "; covariate mode: ", cfg$covariate_mode),
  "",
  "## Top rows",
  "",
  paste(colnames(top), collapse = " | "),
  paste(rep("---", ncol(top)), collapse = " | ")
)
top_rows <- apply(utils::head(top, 40), 1L, function(x) paste(x, collapse = " | "))
writeLines(c(md, top_rows), file.path(cfg$out_dir, "summary.md"))
message("Wrote ", csv)
message("Wrote ", file.path(cfg$out_dir, "summary.md"))
print(utils::head(top, 20), row.names = FALSE)
