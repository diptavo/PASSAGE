# Run conditional PASSAGE H1/H2/H3 on 10x breast Visium Hallmark pathways.
#
# H1 adjusts technical/spatial covariates only.
# H2 additionally adjusts Wu deconvolved cell-type proportions Z_CT.
# H3 additionally adjusts broad background factors V_BG.

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
    seed = 20260523L,
    out_root = "results/passage_10x_hallmark_conditional_wu_perm999",
    deconv_root = "results/wu_deconvolution_10x",
    covariate_mode = "technical_spatial_quadratic",
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    spasset_in_calibration = FALSE,
    datasets = c("Visium_FFPE_Human_Breast_Cancer", "V1_Breast_Cancer_Block_A_Section_1")
  )
  parse_bool <- function(x) {
    x <- tolower(x)
    if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
    if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
    stop("Expected boolean, got: ", x)
  }
  for (arg in args) {
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--cores=", arg)) cfg$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--deconv-root=", arg)) cfg$deconv_root <- sub("^--deconv-root=", "", arg)
    if (grepl("^--covariate-mode=", arg)) cfg$covariate_mode <- sub("^--covariate-mode=", "", arg)
    if (grepl("^--min-pathway-size=", arg)) cfg$min_pathway_size <- as.integer(sub("^--min-pathway-size=", "", arg))
    if (grepl("^--max-pathway-size=", arg)) cfg$max_pathway_size <- as.integer(sub("^--max-pathway-size=", "", arg))
    if (grepl("^--spasset-in-calibration=", arg)) cfg$spasset_in_calibration <- parse_bool(sub("^--spasset-in-calibration=", "", arg))
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
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

prepare_dataset <- function(config, hallmark_pathways, covariate_mode, min_detected_fraction = 0.01, scale_factor = 1e4) {
  counts_raw <- read_10x_h5(file.path(config$root, "filtered_feature_bc_matrix.h5"))
  pos <- read_tissue_positions(file.path(config$root, "spatial"))
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(counts_raw), pos$barcode)
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]
  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)
  counts <- collapse_counts_by_symbol(counts_raw)
  hallmark_genes <- sort(unique(unlist(hallmark_pathways, use.names = FALSE)))
  counts <- counts[intersect(hallmark_genes, rownames(counts)), , drop = FALSE]
  min_spots <- max(20L, ceiling(min_detected_fraction * ncol(counts)))
  counts <- counts[Matrix::rowSums(counts > 0) >= min_spots, , drop = FALSE]
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
  pathways <- pathways[lengths(pathways) >= 5L]
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
    )
  )
}

read_covariate_csv <- function(path, barcodes) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- x$barcode
  mat <- as.matrix(x[, setdiff(names(x), "barcode"), drop = FALSE])
  rownames(mat) <- rn
  mat <- mat[barcodes, , drop = FALSE]
  if (anyNA(mat)) stop("Missing covariate rows in ", path)
  scale(mat)
}

score_one <- function(ii, dat, engine, pre_h1, pre_h2, pre_h3, cfg) {
  pname <- names(dat$pathways)[ii]
  genes <- dat$pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  common_args <- list(
    engine = engine,
    Y = dat$Y,
    pathway = genes,
    weight_schemes = c("equal", "var", "range"),
    gene_names = colnames(dat$Y),
    calibration = "permutation",
    n_perm = cfg$n_perm
  )
  h1 <- do.call(passage_score_test, c(common_args, list(
    precomp = pre_h1,
    seed = cfg$seed + 1000L * ii + 1L,
    run_spasset = cfg$spasset_in_calibration
  )))
  h2 <- do.call(passage_score_test, c(common_args, list(
    precomp = pre_h2,
    seed = cfg$seed + 1000L * ii + 2L,
    run_spasset = FALSE
  )))
  h3 <- do.call(passage_score_test, c(common_args, list(
    precomp = pre_h3,
    seed = cfg$seed + 1000L * ii + 3L,
    run_spasset = FALSE
  )))
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
  decomp <- passage_decomposition(h1, h2, h3)
  pve <- passage_pve(engine, dat$Y, genes, gene_names = colnames(dat$Y))
  data.frame(
    pathway = pname,
    pathway_size = length(genes),
    p_H1 = decomp$p_H1,
    p_H2 = decomp$p_H2,
    p_H3 = decomp$p_H3,
    p_H1_moment = passage_moment_p(h1),
    p_H2_moment = passage_moment_p(h2),
    p_H3_moment = passage_moment_p(h3),
    Q_H1 = decomp$Q_H1,
    Q_H2 = decomp$Q_H2,
    Q_H3 = decomp$Q_H3,
    cell_type_share = decomp$cell_type_share,
    background_share = decomp$background_share,
    pathway_specific_share = decomp$pathway_specific_share,
    R2_cca = pve$summary["R2_cca"],
    PSVS_range = pve$summary["PSVS_range"],
    mean_propSV = pve$summary["mean_propSV"],
    p_spasset = if (!is.null(spasset)) spasset$p else NA_real_,
    spasset_size = if (!is.null(spasset)) spasset$best_size else NA_integer_,
    spasset_genes = if (!is.null(spasset)) paste(spasset$best_genes, collapse = ";") else "",
    elapsed_sec = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
}

run_dataset <- function(dataset_id, hallmark_pathways, cfg) {
  config <- dataset_configs[[dataset_id]]
  message("\n=== Conditional PASSAGE: ", config$label, " ===")
  dat <- prepare_dataset(config, hallmark_pathways, cfg$covariate_mode)
  message("Prepared ", nrow(dat$Y), " spots x ", ncol(dat$Y), " Hallmark genes")
  deconv_dir <- file.path(cfg$deconv_root, dataset_id)
  Z_CT <- read_covariate_csv(file.path(deconv_dir, "cell_type_proportions.csv"), dat$spot_metadata$barcode)
  V_BG <- read_covariate_csv(file.path(deconv_dir, "background_factors.csv"), dat$spot_metadata$barcode)
  message("Loaded Z_CT: ", ncol(Z_CT), " cell types; V_BG: ", ncol(V_BG), " background factors")

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
  X1 <- passage_prepare_design(dat$X, nrow(dat$Y))
  X2 <- cbind(as.matrix(dat$X), as.matrix(Z_CT))
  X3 <- cbind(as.matrix(dat$X), as.matrix(Z_CT), as.matrix(V_BG))
  pre_h1 <- passage_h_precompute(engine, X = X1)
  pre_h2 <- passage_h_precompute(engine, X = passage_prepare_design(X2, nrow(dat$Y)))
  pre_h3 <- passage_h_precompute(engine, X = passage_prepare_design(X3, nrow(dat$Y)))

  keep <- lengths(dat$pathways) >= cfg$min_pathway_size & lengths(dat$pathways) <= cfg$max_pathway_size
  dat$pathways <- dat$pathways[keep]
  message("Scoring ", length(dat$pathways), " pathways for H1/H2/H3 with n_perm=", cfg$n_perm,
          " using ", cfg$cores, " core(s)")
  rows <- parallel::mclapply(
    seq_along(dat$pathways),
    score_one,
    dat = dat,
    engine = engine,
    pre_h1 = pre_h1,
    pre_h2 = pre_h2,
    pre_h3 = pre_h3,
    cfg = cfg,
    mc.cores = cfg$cores,
    mc.preschedule = FALSE
  )
  tbl <- do.call(rbind, rows)
  tbl$fdr_H1 <- stats::p.adjust(tbl$p_H1, method = "BH")
  tbl$fdr_H2 <- stats::p.adjust(tbl$p_H2, method = "BH")
  tbl$fdr_H3 <- stats::p.adjust(tbl$p_H3, method = "BH")
  tbl <- tbl[order(tbl$p_H3, tbl$p_H2, -tbl$pathway_specific_share), , drop = FALSE]
  rownames(tbl) <- NULL

  out_dir <- file.path(cfg$out_root, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- list(
    summary = tbl,
    engine = engine,
    precomp = list(H1 = pre_h1, H2 = pre_h2, H3 = pre_h3),
    Z_CT = Z_CT,
    V_BG = V_BG,
    X_columns = colnames(dat$X),
    Z_CT_columns = colnames(Z_CT),
    V_BG_columns = colnames(V_BG),
    dataset = dataset_id,
    label = config$label,
    n_perm = cfg$n_perm,
    K = cfg$K,
    m = cfg$m,
    covariate_mode = cfg$covariate_mode
  )
  class(fit) <- c("passage_10x_conditional_run", "list")
  saveRDS(fit, file.path(out_dir, "passage_hallmark_conditional_result.rds"))
  saveRDS(dat, file.path(out_dir, "passage_hallmark_conditional_prepared_data.rds"))
  write.csv(tbl, file.path(out_dir, "passage_hallmark_conditional_pathways.csv"), row.names = FALSE)
  write.csv(dat$spot_metadata, file.path(out_dir, "spot_metadata.csv"), row.names = FALSE)

  top <- head(tbl, 15)
  md <- c(
    paste0("# ", config$label, " - Conditional PASSAGE Hallmark"),
    "",
    paste0("- Spots: ", nrow(dat$Y)),
    paste0("- Hallmark genes analyzed: ", ncol(dat$Y)),
    paste0("- Pathways tested: ", nrow(tbl)),
    paste0("- n_perm: ", cfg$n_perm),
    paste0("- Covariate mode: ", cfg$covariate_mode),
    paste0("- H1 covariates: ", paste(colnames(dat$X), collapse = ", ")),
    paste0("- H2 adds cell types: ", paste(colnames(Z_CT), collapse = ", ")),
    paste0("- H3 adds background factors: ", paste(colnames(V_BG), collapse = ", ")),
    "",
    "pathway | p_H1 | p_H2 | p_H3 | fdr_H3 | cell_type_share | background_share | pathway_specific_share | R2_cca | driver_genes",
    "--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---"
  )
  top_rows <- apply(top[, c("pathway", "p_H1", "p_H2", "p_H3", "fdr_H3",
                           "cell_type_share", "background_share", "pathway_specific_share",
                           "R2_cca", "spasset_genes")], 1L, function(r) paste(r, collapse = " | "))
  writeLines(c(md, top_rows), file.path(out_dir, "passage_hallmark_conditional_summary.md"))
  message("Wrote conditional outputs to ", out_dir)
  print(head(tbl[, c("pathway", "p_H1", "p_H2", "p_H3", "fdr_H3",
                    "cell_type_share", "background_share", "pathway_specific_share")], 10), row.names = FALSE)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
bad <- setdiff(cfg$datasets, names(dataset_configs))
if (length(bad) > 0L) stop("Unknown dataset id(s): ", paste(bad, collapse = ", "))
hallmark_pathways <- get_hallmark_pathways()
for (dataset_id in cfg$datasets) {
  run_dataset(dataset_id, hallmark_pathways, cfg)
}
