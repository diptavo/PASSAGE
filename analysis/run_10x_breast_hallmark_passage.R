# Run PASSAGE on two public 10x Visium breast cancer datasets with MSigDB
# Hallmark pathways.
#
# Expected files:
#   data/raw/<dataset>/filtered_feature_bc_matrix.h5
#   data/raw/<dataset>/spatial/tissue_positions_list.csv
#
# Usage:
#   cd /path/to/PASSAGE
#   Rscript scripts/run_10x_breast_hallmark_passage.R
#   Rscript scripts/run_10x_breast_hallmark_passage.R --n-perm=999 --cores=4

suppressPackageStartupMessages({
  library(Matrix)
  library(msigdbr)
  library(rhdf5)
})

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

parse_args <- function(args) {
  cfg <- list(
    n_perm = 999L,
    cores = max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L)),
    K = 6L,
    m = 20L,
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    seed = 20260523L,
    out_root = "results",
    spasset_in_calibration = FALSE,
    covariate_mode = "technical",
    datasets = c("Visium_FFPE_Human_Breast_Cancer", "V1_Breast_Cancer_Block_A_Section_1")
  )
  parse_bool <- function(x) {
    x <- tolower(x)
    if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
    if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
    stop("Expected a boolean value, got: ", x)
  }
  for (arg in args) {
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--cores=", arg)) cfg$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--min-pathway-size=", arg)) cfg$min_pathway_size <- as.integer(sub("^--min-pathway-size=", "", arg))
    if (grepl("^--max-pathway-size=", arg)) cfg$max_pathway_size <- as.integer(sub("^--max-pathway-size=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--spasset-in-calibration=", arg)) {
      cfg$spasset_in_calibration <- parse_bool(sub("^--spasset-in-calibration=", "", arg))
    }
    if (grepl("^--covariate-mode=", arg)) cfg$covariate_mode <- sub("^--covariate-mode=", "", arg)
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
  }
  allowed_modes <- c("technical", "technical_spatial_linear", "technical_spatial_quadratic")
  if (!cfg$covariate_mode %in% allowed_modes) {
    stop("--covariate-mode must be one of: ", paste(allowed_modes, collapse = ", "))
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

get_hallmark_pathways <- function() {
  hall <- msigdbr(species = "Homo sapiens", collection = "H")
  hall$gene_symbol <- toupper(hall$gene_symbol)
  lapply(split(hall$gene_symbol, hall$gs_name), unique)
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
  if (!file.exists(f)) {
    stop("Could not find tissue positions file in ", spatial_dir)
  }
  pos <- read.csv(f, header = header, stringsAsFactors = FALSE)
  if (!header) {
    colnames(pos) <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
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

prepare_dataset <- function(config,
                            hallmark_pathways,
                            covariate_mode = "technical",
                            min_detected_fraction = 0.01,
                            scale_factor = 1e4) {
  h5_path <- file.path(config$root, "filtered_feature_bc_matrix.h5")
  spatial_dir <- file.path(config$root, "spatial")
  if (!file.exists(h5_path)) {
    stop("Missing 10x matrix: ", h5_path)
  }

  counts_raw <- read_10x_h5(h5_path)
  pos <- read_tissue_positions(spatial_dir)
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]

  common <- intersect(colnames(counts_raw), pos$barcode)
  if (length(common) == 0L) {
    stop("No overlapping barcodes between counts and spatial positions")
  }
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]

  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)

  counts <- collapse_counts_by_symbol(counts_raw)
  hallmark_genes <- sort(unique(unlist(hallmark_pathways, use.names = FALSE)))
  keep_hallmark <- intersect(hallmark_genes, rownames(counts))
  counts <- counts[keep_hallmark, , drop = FALSE]

  min_spots <- max(20, ceiling(min_detected_fraction * ncol(counts)))
  detected <- Matrix::rowSums(counts > 0)
  counts <- counts[detected >= min_spots, , drop = FALSE]

  norm <- counts %*% Matrix::Diagonal(x = scale_factor / pmax(lib_size, 1))
  norm@x <- log1p(norm@x)
  Y <- t(as.matrix(norm))

  keep_var <- apply(Y, 2L, stats::var) > 1e-8
  Y <- Y[, keep_var, drop = FALSE]

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

  pathways <- lapply(hallmark_pathways, intersect, y = colnames(Y))
  pathways <- pathways[lengths(pathways) >= 5]

  list(
    Y = Y,
    coords = coords,
    X = X,
    pathways = pathways,
    spot_metadata = data.frame(
      barcode = common,
      lib_size = as.numeric(lib_size),
      n_detected = as.numeric(n_detected),
      pos,
      row.names = NULL
    ),
    n_spots = nrow(Y),
    n_features = ncol(Y),
    n_pathways = length(pathways)
  )
}

score_one_pathway <- function(ii, dat, engine, pre_h1, cfg) {
  pname <- names(dat$pathways)[ii]
  genes <- dat$pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  h1 <- passage_score_test(
    engine = engine,
    Y = dat$Y,
    pathway = genes,
    precomp = pre_h1,
    weight_schemes = c("equal", "var", "range"),
    gene_names = colnames(dat$Y),
    calibration = "permutation",
    n_perm = cfg$n_perm,
    seed = cfg$seed + 1000L * ii,
    run_spasset = cfg$spasset_in_calibration
  )
  h1_driver <- NULL
  if (!cfg$spasset_in_calibration) {
    h1_driver <- passage_score_test(
      engine = engine,
      Y = dat$Y,
      pathway = genes,
      precomp = pre_h1,
      weight_schemes = c("equal", "var", "range"),
      gene_names = colnames(dat$Y),
      calibration = "moment",
      run_burden = FALSE,
      run_spasset = TRUE
    )
  }
  spasset <- if (!is.null(h1$spasset)) h1$spasset else h1_driver$spasset
  pve <- passage_pve(engine, dat$Y, genes, gene_names = colnames(dat$Y))
  elapsed <- proc.time()[["elapsed"]] - t0
  data.frame(
    pathway = pname,
    pathway_size = h1$pathway_size,
    p_H1 = h1$p_omnibus,
    p_H1_moment = h1$p_omnibus_moment,
    p_joint_equal = h1$joint$equal$p_acat,
    Q_H1_equal = h1$joint$equal$Q,
    R2_cca = pve$summary["R2_cca"],
    PSVS_range = pve$summary["PSVS_range"],
    mean_propSV = pve$summary["mean_propSV"],
    p_spasset = if (!is.null(spasset)) spasset$p else NA_real_,
    spasset_size = if (!is.null(spasset)) spasset$best_size else NA_integer_,
    spasset_genes = if (!is.null(spasset)) paste(spasset$best_genes, collapse = ";") else "",
    calibration_components = if (cfg$spasset_in_calibration) {
      "joint+burden+SpASSET"
    } else {
      "joint+burden; SpASSET observed-only"
    },
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

run_passage_dataset <- function(dataset_id, config, hallmark_pathways, cfg) {
  message("\n=== ", config$label, " ===")
  dat <- prepare_dataset(config, hallmark_pathways, covariate_mode = cfg$covariate_mode)
  message("Prepared ", dat$n_spots, " spots, ", dat$n_features,
          " Hallmark genes, ", dat$n_pathways, " pathways")
  message("Covariate mode: ", cfg$covariate_mode, " (", ncol(dat$X), " covariates)")

  range_grid <- passage_default_range_grid(dat$coords, n_grid = 5, min_frac = 0.04, max_frac = 0.50)
  message("Fitting PASSAGE engine")
  engine <- passage_fit_engine_pca(
    Y = dat$Y,
    coords = dat$coords,
    X = dat$X,
    K = cfg$K,
    m = cfg$m,
    range_grid = range_grid,
    kernel = "matern32",
    verbose = TRUE
  )
  pre_h1 <- passage_h_precompute(engine, X = passage_prepare_design(dat$X, nrow(dat$Y)))

  keep <- lengths(dat$pathways) >= cfg$min_pathway_size & lengths(dat$pathways) <= cfg$max_pathway_size
  dat$pathways <- dat$pathways[keep]
  message("Scoring ", length(dat$pathways), " Hallmark pathways with n_perm=", cfg$n_perm,
          " using ", cfg$cores, " core(s)")
  if (!cfg$spasset_in_calibration) {
    message("Permutation calibration uses joint pathway and burden components; SpASSET is computed once per pathway for driver genes.")
  }

  idx <- seq_along(dat$pathways)
  rows <- parallel::mclapply(
    idx,
    score_one_pathway,
    dat = dat,
    engine = engine,
    pre_h1 = pre_h1,
    cfg = cfg,
    mc.cores = cfg$cores,
    mc.preschedule = FALSE
  )
  tbl <- do.call(rbind, rows)
  tbl$fdr_H1 <- stats::p.adjust(tbl$p_H1, method = "BH")
  tbl$fdr_H1_moment <- stats::p.adjust(tbl$p_H1_moment, method = "BH")
  tbl <- tbl[order(tbl$p_H1, tbl$p_H1_moment), , drop = FALSE]
  rownames(tbl) <- NULL

  fit <- list(
    summary = tbl,
    engine = engine,
    precomp = list(H1 = pre_h1),
    ranges = range_grid,
    n_perm = cfg$n_perm,
    K = cfg$K,
    m = cfg$m,
    spasset_in_calibration = cfg$spasset_in_calibration,
    covariate_mode = cfg$covariate_mode,
    X_columns = colnames(dat$X),
    dataset = dataset_id,
    label = config$label
  )
  class(fit) <- c("passage_10x_run", "list")

  out_dir <- file.path(cfg$out_root, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, file.path(out_dir, "passage_hallmark_result.rds"))
  saveRDS(dat, file.path(out_dir, "passage_hallmark_prepared_data.rds"))
  write.csv(tbl, file.path(out_dir, "passage_hallmark_pathways.csv"), row.names = FALSE)
  write.csv(dat$spot_metadata, file.path(out_dir, "spot_metadata.csv"), row.names = FALSE)

  top <- utils::head(tbl, 15)
  md <- c(
    paste0("# ", config$label, " - PASSAGE Hallmark"),
    "",
    paste0("- Spots: ", dat$n_spots),
    paste0("- Hallmark genes analyzed: ", dat$n_features),
    paste0("- Hallmark pathways tested: ", length(dat$pathways)),
    paste0("- K: ", cfg$K),
    paste0("- Vecchia neighbors m: ", cfg$m),
    paste0("- Permutations per pathway: ", cfg$n_perm),
    paste0("- SpASSET included in permutation calibration: ", cfg$spasset_in_calibration),
    paste0("- Covariate mode: ", cfg$covariate_mode),
    paste0("- Covariates adjusted in engine and H1: ", paste(colnames(dat$X), collapse = ", ")),
    "",
    "## Top Hallmark Pathways",
    "",
    paste(c("pathway", "p_H1", "fdr_H1", "p_H1_moment", "pathway_size", "R2_cca", "PSVS_range", "p_spasset", "spasset_size", "spasset_genes"),
          collapse = " | "),
    paste(rep("---", 10), collapse = " | ")
  )
  top_rows <- apply(top[, c("pathway", "p_H1", "fdr_H1", "p_H1_moment", "pathway_size",
                           "R2_cca", "PSVS_range", "p_spasset", "spasset_size", "spasset_genes")], 1L, function(r) {
    paste(r, collapse = " | ")
  })
  writeLines(c(md, top_rows), file.path(out_dir, "passage_hallmark_summary.md"))
  message("Wrote outputs to ", out_dir)
  print(utils::head(tbl[, c("pathway", "p_H1", "fdr_H1", "pathway_size", "R2_cca", "PSVS_range")], 8), row.names = FALSE)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
bad <- setdiff(cfg$datasets, names(dataset_configs))
if (length(bad) > 0L) {
  stop("Unknown dataset id(s): ", paste(bad, collapse = ", "))
}
hallmark_pathways <- get_hallmark_pathways()

for (dataset_id in cfg$datasets) {
  run_passage_dataset(dataset_id, dataset_configs[[dataset_id]], hallmark_pathways, cfg)
}
