# Run conditional PASSAGE H1/H2/H3/H4 on 10x breast Visium MSigDB pathways.
#
# H1 adjusts technical/spatial covariates only.
# H2 additionally adjusts Wu deconvolved cell-type proportions Z_CT.
# H3 additionally adjusts broad background factors V_BG.
# H4 tests whether the residual spatial signature differs across a two-level
# region label. The default region is an epithelial-enrichment median split
# derived from Wu deconvolution, so H4 should be interpreted as exploratory.
#
# The script also computes a per-gene SVG p-value baseline and a pathway-level
# Cauchy/ACAT combined SVG p-value for every tested pathway.

suppressPackageStartupMessages({
  library(Matrix)
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
    engine_rank_method = "fixed",
    engine_method = "pca",
    engine_variance_threshold = 0.90,
    engine_max_K = NA_integer_,
    engine_factor_method = "pca",
    engine_kernel = "matern32",
    engine_n_basis = 60L,
    engine_nmf_iter = 300L,
    engine_alt_iter = 4L,
    engine_smooth_penalty = 1,
    sparse_top_frac = 0.05,
    sparse_min_loadings = 10L,
    seed = 20260523L,
    out_root = "results/passage_10x_kegg_conditional_wu_perm999",
    deconv_root = "results/wu_deconvolution_10x",
    pathway_gmt = file.path("msigdb", "c2.all.v2023.2.Hs.symbols.gmt"),
    pathway_prefix = "KEGG_",
    pathway_label = "KEGG",
    covariate_mode = "technical_spatial_quadratic",
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    max_pathways = NA_integer_,
    min_spot_umi = 500,
    min_spot_genes = 200,
    max_mito_fraction = 0.25,
    max_library_quantile = 0.995,
    h4 = TRUE,
    h4_n_perm = NA_integer_,
    h4_region_method = "epithelial_median",
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
    if (grepl("^--engine-method=", arg)) cfg$engine_method <- sub("^--engine-method=", "", arg)
    if (grepl("^--engine-rank-method=", arg)) cfg$engine_rank_method <- sub("^--engine-rank-method=", "", arg)
    if (grepl("^--engine-variance-threshold=", arg)) cfg$engine_variance_threshold <- as.numeric(sub("^--engine-variance-threshold=", "", arg))
    if (grepl("^--engine-max-K=", arg)) cfg$engine_max_K <- as.integer(sub("^--engine-max-K=", "", arg))
    if (grepl("^--engine-factor-method=", arg)) cfg$engine_factor_method <- sub("^--engine-factor-method=", "", arg)
    if (grepl("^--engine-kernel=", arg)) cfg$engine_kernel <- sub("^--engine-kernel=", "", arg)
    if (grepl("^--engine-n-basis=", arg)) cfg$engine_n_basis <- as.integer(sub("^--engine-n-basis=", "", arg))
    if (grepl("^--engine-nmf-iter=", arg)) cfg$engine_nmf_iter <- as.integer(sub("^--engine-nmf-iter=", "", arg))
    if (grepl("^--engine-alt-iter=", arg)) cfg$engine_alt_iter <- as.integer(sub("^--engine-alt-iter=", "", arg))
    if (grepl("^--engine-smooth-penalty=", arg)) cfg$engine_smooth_penalty <- as.numeric(sub("^--engine-smooth-penalty=", "", arg))
    if (grepl("^--sparse-top-frac=", arg)) cfg$sparse_top_frac <- as.numeric(sub("^--sparse-top-frac=", "", arg))
    if (grepl("^--sparse-min-loadings=", arg)) cfg$sparse_min_loadings <- as.integer(sub("^--sparse-min-loadings=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out-root=", arg)) cfg$out_root <- sub("^--out-root=", "", arg)
    if (grepl("^--deconv-root=", arg)) cfg$deconv_root <- sub("^--deconv-root=", "", arg)
    if (grepl("^--pathway-gmt=", arg)) cfg$pathway_gmt <- sub("^--pathway-gmt=", "", arg)
    if (grepl("^--pathway-prefix=", arg)) cfg$pathway_prefix <- sub("^--pathway-prefix=", "", arg)
    if (grepl("^--pathway-label=", arg)) cfg$pathway_label <- sub("^--pathway-label=", "", arg)
    if (grepl("^--covariate-mode=", arg)) cfg$covariate_mode <- sub("^--covariate-mode=", "", arg)
    if (grepl("^--min-pathway-size=", arg)) cfg$min_pathway_size <- as.integer(sub("^--min-pathway-size=", "", arg))
    if (grepl("^--max-pathway-size=", arg)) cfg$max_pathway_size <- as.integer(sub("^--max-pathway-size=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--min-spot-umi=", arg)) cfg$min_spot_umi <- as.numeric(sub("^--min-spot-umi=", "", arg))
    if (grepl("^--min-spot-genes=", arg)) cfg$min_spot_genes <- as.numeric(sub("^--min-spot-genes=", "", arg))
    if (grepl("^--max-mito-fraction=", arg)) cfg$max_mito_fraction <- as.numeric(sub("^--max-mito-fraction=", "", arg))
    if (grepl("^--max-library-quantile=", arg)) cfg$max_library_quantile <- as.numeric(sub("^--max-library-quantile=", "", arg))
    if (grepl("^--h4=", arg)) cfg$h4 <- parse_bool(sub("^--h4=", "", arg))
    if (grepl("^--h4-n-perm=", arg)) cfg$h4_n_perm <- as.integer(sub("^--h4-n-perm=", "", arg))
    if (grepl("^--h4-region-method=", arg)) cfg$h4_region_method <- sub("^--h4-region-method=", "", arg)
    if (grepl("^--spasset-in-calibration=", arg)) cfg$spasset_in_calibration <- parse_bool(sub("^--spasset-in-calibration=", "", arg))
    if (grepl("^--datasets=", arg)) cfg$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
  }
  if (is.na(cfg$h4_n_perm)) cfg$h4_n_perm <- cfg$n_perm
  if (is.na(cfg$engine_max_K)) cfg$engine_max_K <- cfg$K
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
    genes <- genes[nzchar(genes)]
    out[[nm]] <- genes
  }
  if (length(out) == 0L) {
    stop("No pathways loaded from ", gmt_path, " with prefix '", prefix, "'")
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

prepare_dataset <- function(config, pathways_input, covariate_mode, cfg,
                            min_detected_fraction = 0.01, scale_factor = 1e4) {
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
  max_lib <- if (is.finite(cfg$max_library_quantile) && cfg$max_library_quantile > 0 && cfg$max_library_quantile < 1) {
    as.numeric(stats::quantile(lib_size, probs = cfg$max_library_quantile, na.rm = TRUE, names = FALSE))
  } else Inf
  keep_spots <- lib_size >= cfg$min_spot_umi &
    n_detected >= cfg$min_spot_genes &
    mito_fraction <= cfg$max_mito_fraction &
    lib_size <= max_lib
  qc_summary <- data.frame(
    metric = c(
      "spots_in_tissue_before_qc", "spots_after_qc", "removed_low_umi",
      "removed_low_detected", "removed_high_mito", "removed_high_library",
      "min_spot_umi", "min_spot_genes", "max_mito_fraction",
      "max_library_quantile", "max_library_count"
    ),
    value = c(
      length(common), sum(keep_spots), sum(lib_size < cfg$min_spot_umi),
      sum(n_detected < cfg$min_spot_genes), sum(mito_fraction > cfg$max_mito_fraction),
      sum(lib_size > max_lib), cfg$min_spot_umi, cfg$min_spot_genes,
      cfg$max_mito_fraction, cfg$max_library_quantile, max_lib
    )
  )

  counts_raw <- counts_raw[, keep_spots, drop = FALSE]
  pos <- pos[keep_spots, , drop = FALSE]
  common <- common[keep_spots]
  lib_size <- lib_size[keep_spots]
  n_detected <- n_detected[keep_spots]
  mito_fraction <- mito_fraction[keep_spots]

  counts <- collapse_counts_by_symbol(counts_raw)
  pathway_genes <- sort(unique(unlist(pathways_input, use.names = FALSE)))
  counts <- counts[intersect(pathway_genes, rownames(counts)), , drop = FALSE]
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
  pathways <- lapply(pathways_input, intersect, y = colnames(Y))
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
      mito_fraction = mito_fraction,
      pos,
      row.names = NULL
    ),
    qc_summary = qc_summary
  )
}

read_covariate_csv <- function(path, barcodes, scale_cols = TRUE) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- x$barcode
  mat <- as.matrix(x[, setdiff(names(x), "barcode"), drop = FALSE])
  rownames(mat) <- rn
  mat <- mat[barcodes, , drop = FALSE]
  if (anyNA(mat)) stop("Missing covariate rows in ", path)
  if (scale_cols) scale(mat) else mat
}

make_h4_region <- function(Z_raw, method = "epithelial_median") {
  method <- match.arg(method, c("epithelial_median", "malignant_median", "caf_median"))
  get_col <- function(nm) if (nm %in% colnames(Z_raw)) Z_raw[, nm] else rep(0, nrow(Z_raw))
  score <- switch(
    method,
    epithelial_median = get_col("malignant_epithelial") + get_col("normal_epithelial"),
    malignant_median = get_col("malignant_epithelial"),
    caf_median = get_col("caf")
  )
  cut <- stats::median(score, na.rm = TRUE)
  region <- ifelse(score >= cut, paste0(method, "_high"), paste0(method, "_low"))
  factor(region, levels = c(paste0(method, "_low"), paste0(method, "_high")))
}

passage_h4_region_test <- function(engine, Y, genes, precomp, region,
                                   gene_names, n_perm, seed,
                                   weight_scheme = "equal") {
  P <- passage_resolve_pathway(genes, gene_names)
  R_P <- passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
  A_P <- engine$A[P, , drop = FALSE]
  w <- passage_factor_weights(engine, weight_scheme)
  q_region <- function(region_factor) {
    out <- matrix(0, nrow = engine$K, ncol = 2L)
    colnames(out) <- levels(region_factor)
    for (k in seq_len(engine$K)) {
      vc <- engine$vecchia[[k]]
      region_ord <- region_factor[vc$ord]
      z <- as.numeric(R_P[vc$ord, , drop = FALSE] %*% A_P[, k])
      for (rr in seq_along(levels(region_factor))) {
        mask <- region_ord == levels(region_factor)[rr]
        z_sub <- z[mask]
        K_sub <- engine$K_score[[k]][mask, mask, drop = FALSE]
        out[k, rr] <- as.numeric(crossprod(z_sub, K_sub %*% z_sub)) / max(sum(mask), 1L)
      }
    }
    out
  }
  Q_obs <- q_region(region)
  T_obs <- sum(w * abs(Q_obs[, 2L] - Q_obs[, 1L]))
  if (!is.null(seed)) set.seed(seed)
  T_perm <- numeric(n_perm)
  for (bb in seq_len(n_perm)) {
    region_perm <- factor(sample(region), levels = levels(region))
    Qb <- q_region(region_perm)
    T_perm[bb] <- sum(w * abs(Qb[, 2L] - Qb[, 1L]))
  }
  list(
    p = (1 + sum(T_perm >= T_obs)) / (1 + n_perm),
    T_obs = T_obs,
    Q_region = Q_obs,
    region_levels = levels(region),
    region_sizes = as.integer(table(region)),
    n_perm = n_perm,
    weight_scheme = weight_scheme
  )
}

score_gene_svg <- function(gene, engine, Y, pre_h1, seed) {
  res <- passage_score_test(
    engine = engine,
    Y = Y,
    pathway = gene,
    precomp = pre_h1,
    weight_schemes = c("equal", "var", "range"),
    run_burden = FALSE,
    run_spasset = FALSE,
    gene_names = colnames(Y),
    calibration = "moment",
    seed = seed
  )
  data.frame(
    gene = gene,
    svg_p = res$p_omnibus,
    svg_p_moment = passage_moment_p(res),
    Q_svg_equal = res$joint$equal$Q,
    stringsAsFactors = FALSE
  )
}

compute_gene_svg_table <- function(dat, engine, pre_h1, cfg) {
  genes <- sort(unique(unlist(dat$pathways, use.names = FALSE)))
  message("Running per-gene SVG baseline for ", length(genes), " unique pathway genes using ", cfg$cores, " core(s)")
  rows <- parallel::mclapply(
    seq_along(genes),
    function(ii) score_gene_svg(genes[[ii]], engine, dat$Y, pre_h1, cfg$seed + 500000L + ii),
    mc.cores = cfg$cores,
    mc.preschedule = FALSE
  )
  tbl <- do.call(rbind, rows)
  tbl$svg_fdr <- stats::p.adjust(tbl$svg_p, method = "BH")
  tbl[order(tbl$svg_p, tbl$gene), , drop = FALSE]
}

svg_acat_for_pathway <- function(genes, gene_svg) {
  p <- gene_svg$svg_p[match(genes, gene_svg$gene)]
  p <- p[is.finite(p)]
  if (length(p) == 0L) return(NA_real_)
  passage_acat(p)
}

score_error_row <- function(pname, genes, gene_svg, elapsed_sec, error_message) {
  data.frame(
    pathway = pname,
    pathway_size = length(genes),
    p_H1 = NA_real_,
    p_H2 = NA_real_,
    p_H3 = NA_real_,
    p_H4 = NA_real_,
    svg_acat_p = svg_acat_for_pathway(genes, gene_svg),
    p_H1_moment = NA_real_,
    p_H2_moment = NA_real_,
    p_H3_moment = NA_real_,
    Q_H1 = NA_real_,
    Q_H2 = NA_real_,
    Q_H3 = NA_real_,
    T_H4 = NA_real_,
    cell_type_share = NA_real_,
    background_share = NA_real_,
    pathway_specific_share = NA_real_,
    R2_cca = NA_real_,
    PSVS_range = NA_real_,
    mean_propSV = NA_real_,
    p_spasset = NA_real_,
    spasset_size = NA_integer_,
    spasset_genes = "",
    elapsed_sec = elapsed_sec,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

score_one_impl <- function(ii, dat, engine, pre_h1, pre_h2, pre_h3, pre_h4, region_h4, gene_svg, cfg) {
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
  h4 <- NULL
  if (isTRUE(cfg$h4)) {
    h4 <- passage_region_enrichment_test(
      engine = engine,
      Y = dat$Y,
      pathway = genes,
      precomp = pre_h4,
      region = region_h4,
      gene_names = colnames(dat$Y),
      n_perm = cfg$h4_n_perm,
      seed = cfg$seed + 1000L * ii + 4L
    )
  }
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
  row <- data.frame(
    pathway = pname,
    pathway_size = length(genes),
    p_H1 = decomp$p_H1,
    p_H2 = decomp$p_H2,
    p_H3 = decomp$p_H3,
    p_H4 = if (!is.null(h4)) h4$p else NA_real_,
    svg_acat_p = svg_acat_for_pathway(genes, gene_svg),
    p_H1_moment = passage_moment_p(h1),
    p_H2_moment = passage_moment_p(h2),
    p_H3_moment = passage_moment_p(h3),
    Q_H1 = decomp$Q_H1,
    Q_H2 = decomp$Q_H2,
    Q_H3 = decomp$Q_H3,
    T_H4 = if (!is.null(h4)) h4$T_obs else NA_real_,
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
    error_message = "",
    stringsAsFactors = FALSE
  )
  msg <- sprintf(
    "PATHWAY_DONE\t%s\tp_H1=%.4g\tp_H2=%.4g\tp_H3=%.4g\tp_H4=%.4g\tsvg_acat=%.4g\tspecific=%.3f\telapsed=%.1fs",
    pname, row$p_H1, row$p_H2, row$p_H3, row$p_H4, row$svg_acat_p,
    row$pathway_specific_share, row$elapsed_sec
  )
  message(msg)
  row
}

score_one <- function(ii, dat, engine, pre_h1, pre_h2, pre_h3, pre_h4, region_h4, gene_svg, cfg) {
  pname <- names(dat$pathways)[ii]
  genes <- dat$pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  tryCatch(
    score_one_impl(ii, dat, engine, pre_h1, pre_h2, pre_h3, pre_h4, region_h4, gene_svg, cfg),
    error = function(e) {
      elapsed_sec <- proc.time()[["elapsed"]] - t0
      msg <- conditionMessage(e)
      message(sprintf("PATHWAY_ERROR\t%s\t%s\telapsed=%.1fs", pname, msg, elapsed_sec))
      score_error_row(pname, genes, gene_svg, elapsed_sec, msg)
    }
  )
}

run_dataset <- function(dataset_id, pathways_input, cfg) {
  config <- dataset_configs[[dataset_id]]
  message("\n=== Conditional PASSAGE ", cfg$pathway_label, ": ", config$label, " ===")
  dat <- prepare_dataset(config, pathways_input, cfg$covariate_mode, cfg)
  message("Prepared ", nrow(dat$Y), " spots x ", ncol(dat$Y), " ", cfg$pathway_label, " genes after QC")
  print(dat$qc_summary, row.names = FALSE)

  deconv_dir <- file.path(cfg$deconv_root, dataset_id)
  Z_raw <- read_covariate_csv(file.path(deconv_dir, "cell_type_proportions.csv"), dat$spot_metadata$barcode, scale_cols = FALSE)
  Z_CT <- read_covariate_csv(file.path(deconv_dir, "cell_type_proportions.csv"), dat$spot_metadata$barcode)
  V_BG <- read_covariate_csv(file.path(deconv_dir, "background_factors.csv"), dat$spot_metadata$barcode)
  message("Loaded Z_CT: ", ncol(Z_CT), " cell types; V_BG: ", ncol(V_BG), " background factors")
  region_h4 <- make_h4_region(Z_raw, cfg$h4_region_method)
  message("H4 region: ", paste(names(table(region_h4)), as.integer(table(region_h4)), sep = "=", collapse = "; "))

  range_grid <- passage_default_range_grid(dat$coords, n_grid = 5, min_frac = 0.04, max_frac = 0.50)
  message("Fitting PASSAGE engine: ", cfg$engine_method)
  engine <- switch(
    cfg$engine_method,
    pca = passage_fit_engine_pca(
      Y = dat$Y,
      coords = dat$coords,
      X = dat$X,
      K = cfg$K,
      rank_method = cfg$engine_rank_method,
      variance_threshold = cfg$engine_variance_threshold,
      max_K = cfg$engine_max_K,
      factor_method = cfg$engine_factor_method,
      sparse_top_frac = cfg$sparse_top_frac,
      sparse_min_loadings = cfg$sparse_min_loadings,
      m = cfg$m,
      range_grid = range_grid,
      kernel = cfg$engine_kernel,
      verbose = TRUE
    ),
    spatial_basis = passage_fit_engine_spatial_basis(
      Y = dat$Y,
      coords = dat$coords,
      X = dat$X,
      K = cfg$K,
      rank_method = cfg$engine_rank_method,
      variance_threshold = cfg$engine_variance_threshold,
      max_K = cfg$engine_max_K,
      n_basis = cfg$engine_n_basis,
      m = cfg$m,
      range_grid = range_grid,
      kernel = cfg$engine_kernel,
      verbose = TRUE
    ),
    smoothed_pca = passage_fit_engine_smoothed_pca(
      Y = dat$Y,
      coords = dat$coords,
      X = dat$X,
      K = cfg$K,
      rank_method = cfg$engine_rank_method,
      variance_threshold = cfg$engine_variance_threshold,
      max_K = cfg$engine_max_K,
      n_basis = cfg$engine_n_basis,
      m = cfg$m,
      range_grid = range_grid,
      kernel = cfg$engine_kernel,
      verbose = TRUE
    ),
    nmf = passage_fit_engine_nmf(
      Y = dat$Y,
      coords = dat$coords,
      X = dat$X,
      K = cfg$K,
      n_iter = cfg$engine_nmf_iter,
      n_basis = cfg$engine_n_basis,
      m = cfg$m,
      range_grid = range_grid,
      kernel = cfg$engine_kernel,
      verbose = TRUE
    ),
    alternating_gp = passage_fit_engine_alternating_gp(
      Y = dat$Y,
      coords = dat$coords,
      X = dat$X,
      K = cfg$K,
      rank_method = cfg$engine_rank_method,
      variance_threshold = cfg$engine_variance_threshold,
      max_K = cfg$engine_max_K,
      n_iter = cfg$engine_alt_iter,
      smooth_penalty = cfg$engine_smooth_penalty,
      m = cfg$m,
      range_grid = range_grid,
      kernel = cfg$engine_kernel,
      verbose = TRUE
    ),
    stop("Unsupported --engine-method: ", cfg$engine_method)
  )

  X1 <- passage_prepare_design(dat$X, nrow(dat$Y))
  X2 <- cbind(as.matrix(dat$X), as.matrix(Z_CT))
  X3 <- cbind(as.matrix(dat$X), as.matrix(Z_CT), as.matrix(V_BG))
  X4 <- cbind(as.matrix(dat$X), as.matrix(Z_CT), h4_region_high = as.numeric(region_h4 == levels(region_h4)[2L]))
  pre_h1 <- passage_h_precompute(engine, X = X1)
  pre_h2 <- passage_h_precompute(engine, X = passage_prepare_design(X2, nrow(dat$Y)))
  pre_h3 <- passage_h_precompute(engine, X = passage_prepare_design(X3, nrow(dat$Y)))
  pre_h4 <- passage_h_precompute(engine, X = passage_prepare_design(X4, nrow(dat$Y)))

  keep <- lengths(dat$pathways) >= cfg$min_pathway_size & lengths(dat$pathways) <= cfg$max_pathway_size
  dat$pathways <- dat$pathways[keep]
  if (is.finite(cfg$max_pathways)) {
    dat$pathways <- dat$pathways[seq_len(min(length(dat$pathways), cfg$max_pathways))]
  }
  message("Retained ", length(dat$pathways), " ", cfg$pathway_label, " pathways after size filter")

  out_dir <- file.path(cfg$out_root, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gene_svg_path <- file.path(out_dir, "gene_svg_pvalues.csv")
  if (file.exists(gene_svg_path)) {
    message("Reusing existing per-gene SVG baseline: ", gene_svg_path)
    gene_svg <- read.csv(gene_svg_path, stringsAsFactors = FALSE)
  } else {
    gene_svg <- compute_gene_svg_table(dat, engine, pre_h1, cfg)
    write.csv(gene_svg, gene_svg_path, row.names = FALSE)
  }

  message("Scoring ", length(dat$pathways), " pathways for H1/H2/H3",
          if (cfg$h4) "/H4" else "", " with n_perm=", cfg$n_perm,
          if (cfg$h4) paste0(" and h4_n_perm=", cfg$h4_n_perm) else "",
          " using ", cfg$cores, " core(s)")
  rows <- parallel::mclapply(
    seq_along(dat$pathways),
    score_one,
    dat = dat,
    engine = engine,
    pre_h1 = pre_h1,
    pre_h2 = pre_h2,
    pre_h3 = pre_h3,
    pre_h4 = pre_h4,
    region_h4 = region_h4,
    gene_svg = gene_svg,
    cfg = cfg,
    mc.cores = cfg$cores,
    mc.preschedule = FALSE
  )
  tbl <- do.call(rbind, rows)
  tbl$fdr_H1 <- stats::p.adjust(tbl$p_H1, method = "BH")
  tbl$fdr_H2 <- stats::p.adjust(tbl$p_H2, method = "BH")
  tbl$fdr_H3 <- stats::p.adjust(tbl$p_H3, method = "BH")
  tbl$fdr_H4 <- stats::p.adjust(tbl$p_H4, method = "BH")
  tbl$svg_acat_fdr <- stats::p.adjust(tbl$svg_acat_p, method = "BH")
  tbl <- tbl[order(tbl$p_H3, tbl$p_H2, tbl$p_H4, -tbl$pathway_specific_share,
                   na.last = TRUE), , drop = FALSE]
  rownames(tbl) <- NULL

  fit <- list(
    summary = tbl,
    engine = engine,
    precomp = list(H1 = pre_h1, H2 = pre_h2, H3 = pre_h3, H4 = pre_h4),
    Z_CT = Z_CT,
    V_BG = V_BG,
    H4_region = region_h4,
    gene_svg = gene_svg,
    X_columns = colnames(dat$X),
    Z_CT_columns = colnames(Z_CT),
    V_BG_columns = colnames(V_BG),
    dataset = dataset_id,
    label = config$label,
    n_perm = cfg$n_perm,
    h4_n_perm = cfg$h4_n_perm,
    K = cfg$K,
    engine_K = engine$K,
    engine_method = cfg$engine_method,
    engine_rank_method = cfg$engine_rank_method,
    engine_variance_threshold = cfg$engine_variance_threshold,
    engine_factor_method = cfg$engine_factor_method,
    engine_kernel = cfg$engine_kernel,
    engine_rank_info = engine$rank_info,
    m = cfg$m,
    covariate_mode = cfg$covariate_mode
  )
  class(fit) <- c("passage_10x_conditional_run", "list")
  saveRDS(fit, file.path(out_dir, "passage_msigdb_conditional_result.rds"))
  saveRDS(dat, file.path(out_dir, "passage_msigdb_conditional_prepared_data.rds"))
  write.csv(tbl, file.path(out_dir, "passage_msigdb_conditional_pathways.csv"), row.names = FALSE)
  write.csv(dat$spot_metadata, file.path(out_dir, "spot_metadata.csv"), row.names = FALSE)
  write.csv(dat$qc_summary, file.path(out_dir, "qc_summary.csv"), row.names = FALSE)
  write.csv(data.frame(barcode = dat$spot_metadata$barcode, H4_region = as.character(region_h4)),
            file.path(out_dir, "h4_region_labels.csv"), row.names = FALSE)

  top <- head(tbl, 15)
  md <- c(
    paste0("# ", config$label, " - Conditional PASSAGE ", cfg$pathway_label),
    "",
    paste0("- Spots: ", nrow(dat$Y)),
    paste0("- ", cfg$pathway_label, " genes analyzed: ", ncol(dat$Y)),
    paste0("- Pathways tested: ", nrow(tbl)),
    paste0("- n_perm: ", cfg$n_perm),
    paste0("- h4_n_perm: ", cfg$h4_n_perm),
    paste0("- Engine: ", cfg$engine_method,
           if (cfg$engine_method == "pca") paste0(":", cfg$engine_factor_method) else "",
           " / ", cfg$engine_rank_method,
           " / ", cfg$engine_kernel, " / K=", engine$K,
           " / cumulative residual variance=", round(engine$rank_info$cumulative_variance, 4)),
    paste0("- Covariate mode: ", cfg$covariate_mode),
    paste0("- H1 covariates: ", paste(colnames(dat$X), collapse = ", ")),
    paste0("- H2 adds cell types: ", paste(colnames(Z_CT), collapse = ", ")),
    paste0("- H3 adds background factors: ", paste(colnames(V_BG), collapse = ", ")),
    paste0("- H4 region: ", cfg$h4_region_method, " (", paste(names(table(region_h4)), as.integer(table(region_h4)), sep = "=", collapse = "; "), ")"),
    "",
    "pathway | p_H1 | p_H2 | p_H3 | p_H4 | svg_acat_p | fdr_H3 | fdr_H4 | svg_acat_fdr | cell_type_share | background_share | pathway_specific_share | R2_cca | driver_genes",
    "--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---"
  )
  top_rows <- apply(top[, c("pathway", "p_H1", "p_H2", "p_H3", "p_H4", "svg_acat_p",
                           "fdr_H3", "fdr_H4", "svg_acat_fdr",
                           "cell_type_share", "background_share", "pathway_specific_share",
                           "R2_cca", "spasset_genes")], 1L, function(r) paste(r, collapse = " | "))
  writeLines(c(md, top_rows), file.path(out_dir, "passage_msigdb_conditional_summary.md"))
  message("Wrote conditional outputs to ", out_dir)
  print(head(tbl[, c("pathway", "p_H1", "p_H2", "p_H3", "p_H4", "svg_acat_p",
                    "fdr_H3", "fdr_H4", "svg_acat_fdr",
                    "cell_type_share", "background_share", "pathway_specific_share")], 10), row.names = FALSE)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
bad <- setdiff(cfg$datasets, names(dataset_configs))
if (length(bad) > 0L) stop("Unknown dataset id(s): ", paste(bad, collapse = ", "))
pathways_input <- read_gmt_pathways(cfg$pathway_gmt, cfg$pathway_prefix)
message("Loaded ", length(pathways_input), " ", cfg$pathway_label, " pathways from ", cfg$pathway_gmt)
for (dataset_id in cfg$datasets) {
  run_dataset(dataset_id, pathways_input, cfg)
}
