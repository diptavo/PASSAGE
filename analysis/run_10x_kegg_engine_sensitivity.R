# Engine sensitivity analysis for the KEGG pathways already streamed in run.log.
#
# This is intentionally exploratory. It compares rank selection, low-rank factor
# extraction, and Matérn kernels on a fixed set of already-completed pathways.

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
})

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

parse_args <- function(args) {
  cfg <- list(
    dataset = "Visium_FFPE_Human_Breast_Cancer",
    run_log = file.path("results", "passage_10x_kegg_conditional_wu_perm999", "run.log"),
    out_root = file.path("results", "passage_10x_kegg_engine_sensitivity"),
    deconv_root = "results/wu_deconvolution_10x",
    pathway_gmt = file.path("msigdb", "c2.all.v2023.2.Hs.symbols.gmt"),
    pathway_prefix = "KEGG_",
    n_perm = 99L,
    cores = max(1L, min(2L, parallel::detectCores(logical = FALSE) - 1L)),
    K = 6L,
    auto_max_K = 12L,
    auto_variance = 0.90,
    m = 20L,
    seed = 20260524L,
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    min_spot_umi = 500,
    min_spot_genes = 200,
    max_mito_fraction = 0.25,
    max_library_quantile = 0.995
  )
  for (arg in args) {
    if (grepl("^--dataset=", arg)) cfg$dataset <- sub("^--dataset=", "", arg)
    if (grepl("^--run-log=", arg)) cfg$run_log <- sub("^--run-log=", "", arg)
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--deconv-root=", arg)) cfg$deconv_root <- sub("^--deconv-root=", "", arg)
    if (grepl("^--pathway-gmt=", arg)) cfg$pathway_gmt <- sub("^--pathway-gmt=", "", arg)
    if (grepl("^--pathway-prefix=", arg)) cfg$pathway_prefix <- sub("^--pathway-prefix=", "", arg)
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--cores=", arg)) cfg$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--auto-max-K=", arg)) cfg$auto_max_K <- as.integer(sub("^--auto-max-K=", "", arg))
    if (grepl("^--auto-variance=", arg)) cfg$auto_variance <- as.numeric(sub("^--auto-variance=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
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

read_completed_pathways <- function(run_log) {
  if (!file.exists(run_log)) stop("run_log not found: ", run_log)
  lines <- readLines(run_log, warn = FALSE)
  lines <- lines[startsWith(lines, "PATHWAY_DONE\t")]
  unique(vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 2L))
}

read_gmt_pathways <- function(gmt_path, prefix = NULL) {
  if (!file.exists(gmt_path)) stop("GMT file not found: ", gmt_path)
  lines <- readLines(gmt_path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) next
    nm <- fields[[1L]]
    if (!is.null(prefix) && nzchar(prefix) && !startsWith(nm, prefix)) next
    genes <- toupper(unique(fields[-c(1L, 2L)]))
    out[[nm]] <- genes[nzchar(genes)]
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

prepare_dataset <- function(config, pathways_input, cfg, scale_factor = 1e4) {
  counts_raw <- read_10x_h5(file.path(config$root, "filtered_feature_bc_matrix.h5"))
  pos <- read_tissue_positions(file.path(config$root, "spatial"))
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  common <- intersect(colnames(counts_raw), pos$barcode)
  counts_raw <- counts_raw[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]

  lib_size <- Matrix::colSums(counts_raw)
  n_detected <- Matrix::colSums(counts_raw > 0)
  mito_rows <- grepl("^MT-", toupper(rownames(counts_raw)))
  mito_counts <- if (any(mito_rows)) Matrix::colSums(counts_raw[mito_rows, , drop = FALSE]) else rep(0, ncol(counts_raw))
  mito_fraction <- as.numeric(mito_counts / pmax(lib_size, 1))
  max_lib <- as.numeric(stats::quantile(lib_size, probs = cfg$max_library_quantile, na.rm = TRUE, names = FALSE))
  keep <- lib_size >= cfg$min_spot_umi &
    n_detected >= cfg$min_spot_genes &
    mito_fraction <= cfg$max_mito_fraction &
    lib_size <= max_lib

  counts_raw <- counts_raw[, keep, drop = FALSE]
  pos <- pos[keep, , drop = FALSE]
  common <- common[keep]
  lib_size <- lib_size[keep]
  n_detected <- n_detected[keep]
  mito_fraction <- mito_fraction[keep]

  counts <- collapse_counts_by_symbol(counts_raw)
  pathway_genes <- sort(unique(unlist(pathways_input, use.names = FALSE)))
  counts <- counts[intersect(pathway_genes, rownames(counts)), , drop = FALSE]
  min_spots <- max(20L, ceiling(0.01 * ncol(counts)))
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
  sx <- as.numeric(scale(coords[, 1L]))
  sy <- as.numeric(scale(coords[, 2L]))
  X <- cbind(
    log_umi = as.numeric(scale(log1p(lib_size))),
    log_detected = as.numeric(scale(log1p(n_detected))),
    spatial_x = sx,
    spatial_y = sy,
    spatial_x2 = as.numeric(scale(sx^2)),
    spatial_y2 = as.numeric(scale(sy^2)),
    spatial_xy = as.numeric(scale(sx * sy))
  )
  rownames(X) <- common
  pathways <- lapply(pathways_input, intersect, y = colnames(Y))
  list(
    Y = Y,
    coords = coords,
    X = X,
    pathways = pathways,
    spot_metadata = data.frame(
      barcode = common,
      lib_size = as.numeric(lib_size),
      n_detected = as.numeric(n_detected),
      mito_fraction = mito_fraction,
      pos,
      row.names = NULL
    )
  )
}

read_covariate_csv <- function(path, barcodes) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(x[, setdiff(names(x), "barcode"), drop = FALSE])
  rownames(mat) <- x$barcode
  mat <- mat[barcodes, , drop = FALSE]
  if (anyNA(mat)) stop("Missing covariate rows in ", path)
  scale(mat)
}

score_pathway <- function(ii, pathways, dat, engine, pre_h1, pre_h2, pre_h3, cfg, variant) {
  pname <- names(pathways)[ii]
  genes <- pathways[[ii]]
  common_args <- list(
    engine = engine,
    Y = dat$Y,
    pathway = genes,
    weight_schemes = c("equal", "var", "range"),
    gene_names = colnames(dat$Y),
    calibration = "permutation",
    n_perm = cfg$n_perm,
    run_spasset = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  h1 <- do.call(passage_score_test, c(common_args, list(precomp = pre_h1, seed = cfg$seed + ii * 1000L + 1L)))
  h2 <- do.call(passage_score_test, c(common_args, list(precomp = pre_h2, seed = cfg$seed + ii * 1000L + 2L)))
  h3 <- do.call(passage_score_test, c(common_args, list(precomp = pre_h3, seed = cfg$seed + ii * 1000L + 3L)))
  d <- passage_decomposition(h1, h2, h3)
  row <- data.frame(
    variant = variant$name,
    pathway = pname,
    pathway_size = length(genes),
    p_H1 = d$p_H1,
    p_H2 = d$p_H2,
    p_H3 = d$p_H3,
    Q_H1 = d$Q_H1,
    Q_H2 = d$Q_H2,
    Q_H3 = d$Q_H3,
    cell_type_share = d$cell_type_share,
    background_share = d$background_share,
    pathway_specific_share = d$pathway_specific_share,
    elapsed_sec = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "SENS_DONE\t%s\t%s\tp_H1=%.4g\tp_H2=%.4g\tp_H3=%.4g\tspecific=%.3f\telapsed=%.1fs",
    variant$name, pname, row$p_H1, row$p_H2, row$p_H3, row$pathway_specific_share, row$elapsed_sec
  ))
  row
}

run_variant <- function(variant, dat, Z_CT, V_BG, pathways, cfg, out_dir) {
  message("\n=== ENGINE_VARIANT ", variant$name, " ===")
  range_grid <- passage_default_range_grid(dat$coords, n_grid = 5, min_frac = 0.04, max_frac = 0.50)
  engine <- passage_fit_engine_pca(
    Y = dat$Y,
    coords = dat$coords,
    X = dat$X,
    K = variant$K,
    rank_method = variant$rank_method,
    variance_threshold = cfg$auto_variance,
    max_K = variant$max_K,
    min_K = 2L,
    factor_method = variant$factor_method,
    m = cfg$m,
    range_grid = range_grid,
    kernel = variant$kernel,
    verbose = TRUE
  )
  X1 <- passage_prepare_design(dat$X, nrow(dat$Y))
  X2 <- passage_prepare_design(cbind(as.matrix(dat$X), as.matrix(Z_CT)), nrow(dat$Y))
  X3 <- passage_prepare_design(cbind(as.matrix(dat$X), as.matrix(Z_CT), as.matrix(V_BG)), nrow(dat$Y))
  pre_h1 <- passage_h_precompute(engine, X = X1)
  pre_h2 <- passage_h_precompute(engine, X = X2)
  pre_h3 <- passage_h_precompute(engine, X = X3)
  rows <- parallel::mclapply(
    seq_along(pathways),
    score_pathway,
    pathways = pathways,
    dat = dat,
    engine = engine,
    pre_h1 = pre_h1,
    pre_h2 = pre_h2,
    pre_h3 = pre_h3,
    cfg = cfg,
    variant = variant,
    mc.cores = cfg$cores,
    mc.preschedule = FALSE
  )
  tbl <- do.call(rbind, rows)
  tbl$engine_K <- engine$K
  tbl$rank_method <- variant$rank_method
  tbl$factor_method <- engine$factor_method
  tbl$kernel <- engine$kernel
  tbl$rank_cumulative_variance <- engine$rank_info$cumulative_variance
  tbl$rank_threshold_reached <- engine$rank_info$threshold_reached
  tbl$fdr_H1 <- stats::p.adjust(tbl$p_H1, method = "BH")
  tbl$fdr_H2 <- stats::p.adjust(tbl$p_H2, method = "BH")
  tbl$fdr_H3 <- stats::p.adjust(tbl$p_H3, method = "BH")
  write.csv(tbl, file.path(out_dir, paste0("variant_", variant$name, ".csv")), row.names = FALSE)
  saveRDS(
    list(summary = tbl, engine_rank_info = engine$rank_info, theta = engine$theta),
    file.path(out_dir, paste0("variant_", variant$name, ".rds"))
  )
  message(sprintf(
    "VARIANT_DONE\t%s\tK=%d\tcumvar=%.3f\tthreshold_reached=%s\tH3_FDR05=%d",
    variant$name, engine$K, engine$rank_info$cumulative_variance,
    engine$rank_info$threshold_reached, sum(tbl$fdr_H3 < 0.05, na.rm = TRUE)
  ))
  tbl
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
if (!cfg$dataset %in% names(dataset_configs)) stop("Unknown dataset: ", cfg$dataset)
completed <- read_completed_pathways(cfg$run_log)
pathways_all <- read_gmt_pathways(cfg$pathway_gmt, cfg$pathway_prefix)
pathways_all <- pathways_all[intersect(completed, names(pathways_all))]
if (length(pathways_all) == 0L) stop("No completed pathways found in GMT")

config <- dataset_configs[[cfg$dataset]]
message("Preparing ", config$label)
dat <- prepare_dataset(config, pathways_all, cfg)
pathways <- lapply(pathways_all, intersect, y = colnames(dat$Y))
keep <- lengths(pathways) >= cfg$min_pathway_size & lengths(pathways) <= cfg$max_pathway_size
pathways <- pathways[keep]
message("Sensitivity pathway set: ", length(pathways), " pathways")

deconv_dir <- file.path(cfg$deconv_root, cfg$dataset)
Z_CT <- read_covariate_csv(file.path(deconv_dir, "cell_type_proportions.csv"), dat$spot_metadata$barcode)
V_BG <- read_covariate_csv(file.path(deconv_dir, "background_factors.csv"), dat$spot_metadata$barcode)

variants <- list(
  list(name = "fixed6_pca_matern32", K = cfg$K, max_K = cfg$K, rank_method = "fixed", factor_method = "pca", kernel = "matern32"),
  list(name = "auto90_pca_matern32", K = cfg$K, max_K = cfg$auto_max_K, rank_method = "variance", factor_method = "pca", kernel = "matern32"),
  list(name = "auto90_varimax_matern32", K = cfg$K, max_K = cfg$auto_max_K, rank_method = "variance", factor_method = "varimax", kernel = "matern32"),
  list(name = "auto90_nmf_matern32", K = cfg$K, max_K = cfg$auto_max_K, rank_method = "variance", factor_method = "nmf", kernel = "matern32"),
  list(name = "auto90_pca_matern12", K = cfg$K, max_K = cfg$auto_max_K, rank_method = "variance", factor_method = "pca", kernel = "matern12"),
  list(name = "auto90_pca_matern52", K = cfg$K, max_K = cfg$auto_max_K, rank_method = "variance", factor_method = "pca", kernel = "matern52")
)

out_dir <- file.path(cfg$out_root, cfg$dataset)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(completed, file.path(out_dir, "completed_pathways_from_run_log.txt"))
all_results <- lapply(variants, run_variant, dat = dat, Z_CT = Z_CT, V_BG = V_BG,
                      pathways = pathways, cfg = cfg, out_dir = out_dir)
combined <- do.call(rbind, all_results)
write.csv(combined, file.path(out_dir, "engine_sensitivity_all_variants.csv"), row.names = FALSE)

summary <- aggregate(
  cbind(H1_sig = combined$fdr_H1 < 0.05, H2_sig = combined$fdr_H2 < 0.05,
        H3_sig = combined$fdr_H3 < 0.05, rank_cumulative_variance = combined$rank_cumulative_variance) ~
    variant + engine_K + rank_method + factor_method + kernel + rank_threshold_reached,
  data = combined,
  FUN = function(x) if (is.logical(x)) sum(x, na.rm = TRUE) else unique(round(x, 4))[1]
)
write.csv(summary, file.path(out_dir, "engine_sensitivity_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
