# Real-data-derived residual-spatial null calibration for PASSAGE.
#
# For each prepared breast Visium dataset, this script residualizes expression
# on the prepared covariate matrix, permutes residual rows with one common spot
# permutation across all genes, adds the fitted covariate mean back, and then
# runs PASSAGE. This preserves gene covariance and covariate effects while
# breaking residual spatial alignment with tissue coordinates.
#
# Usage:
#   Rscript scripts/calibrate_passage_real_residual_null.R \
#     --n-reps=5 --n-perm=19 --max-pathways=4 \
#     --out-dir=results/passage_real_residual_null_calibration_quick

parse_args <- function(args) {
  cfg <- list(
    seed = 20260730L,
    n_reps = 5L,
    n_perm = 19L,
    max_pathways = 4L,
    K = 3L,
    m = 8L,
    n_folds = 4L,
    out_dir = file.path("results", "passage_real_residual_null_calibration_quick"),
    modes = c("none", "spot_crossfit", "pathway_holdout")
  )
  for (arg in args) {
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--n-reps=", arg)) cfg$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--modes=", arg)) cfg$modes <- strsplit(sub("^--modes=", "", arg), ",", fixed = TRUE)[[1]]
  }
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

prepared_paths <- function() {
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_perm999_fastcal/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_prepared_data.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_perm999_fastcal/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_prepared_data.rds"
  )
}

select_pathways <- function(pathways, max_pathways) {
  preferred <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA"
  )
  keep <- c(intersect(preferred, names(pathways)), setdiff(names(pathways), preferred))
  pathways[keep[seq_len(min(max_pathways, length(keep)))]]
}

make_residual_null_y <- function(Y, X, seed) {
  set.seed(seed)
  Xd <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  qrx <- qr(Xd)
  R <- passage_residualize_with_qr(Y, qrx)
  Yhat <- Y - R
  Yhat + R[sample.int(nrow(R)), , drop = FALSE]
}

run_mode <- function(Y, coords, X, pathways, mode, rep_id, cfg) {
  range_grid <- passage_default_range_grid(coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
  fit_args <- list(K = cfg$K, m = cfg$m, range_grid = range_grid,
                   kernel = "matern32", ordering = "none", verbose = FALSE)
  t0 <- proc.time()[["elapsed"]]
  fit <- if (mode == "none") {
    passage_run(Y, coords, pathways, X = X, K = cfg$K, m = cfg$m,
                range_grid = range_grid, hypotheses = "H1",
                calibration = "permutation", n_perm = cfg$n_perm,
                seed = cfg$seed + rep_id, verbose = FALSE)
  } else if (mode == "spot_crossfit") {
    passage_run(Y, coords, pathways, X = X, K = cfg$K, m = cfg$m,
                range_grid = range_grid, hypotheses = "H1",
                calibration = "permutation", n_perm = cfg$n_perm,
                seed = cfg$seed + rep_id, verbose = FALSE,
                anti_dip = "spot_crossfit",
                anti_dip_args = list(n_folds = cfg$n_folds, seed = cfg$seed + rep_id, fit_args = fit_args))
  } else if (mode == "pathway_holdout") {
    passage_run(Y, coords, pathways, X = X, K = cfg$K, m = cfg$m,
                range_grid = range_grid, hypotheses = "H1",
                calibration = "permutation", n_perm = cfg$n_perm,
                seed = cfg$seed + rep_id, verbose = FALSE,
                anti_dip = "pathway_holdout",
                anti_dip_args = list(fit_args = fit_args))
  } else {
    stop("unknown mode: ", mode)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  tbl <- fit$summary
  tbl$mode <- mode
  tbl$replicate <- rep_id
  tbl$elapsed_sec <- elapsed
  tbl[, c("replicate", "mode", "pathway", "pathway_size", "p_H1", "p_H1_moment",
          "R2_cca", "PSVS_range", "mean_propSV", "elapsed_sec")]
}

summarize_results <- function(x) {
  x$p <- x$p_H1
  rows <- list()
  ii <- 0L
  for (dataset in unique(x$dataset)) {
    for (mode in unique(x$mode)) {
      z <- x[x$dataset == dataset & x$mode == mode, , drop = FALSE]
      for (alpha in c(0.10, 0.05, 0.01)) {
        p <- z$p[is.finite(z$p)]
        rate <- mean(p <= alpha)
        se <- sqrt(rate * (1 - rate) / max(1L, length(p)))
        ii <- ii + 1L
        rows[[ii]] <- data.frame(
          dataset = dataset,
          mode = mode,
          alpha = alpha,
          n = length(p),
          reject_rate = rate,
          mc_se = se,
          ci95_low = max(0, rate - 1.96 * se),
          ci95_high = min(1, rate + 1.96 * se),
          median_p = stats::median(p),
          min_p = min(p),
          median_elapsed_sec = stats::median(z$elapsed_sec),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
ii <- 0L
for (dataset in names(prepared_paths())) {
  message("loading ", dataset)
  dat <- readRDS(prepared_paths()[[dataset]])
  pathways <- select_pathways(dat$pathways, cfg$max_pathways)
  for (rr in seq_len(cfg$n_reps)) {
    message("dataset=", dataset, " residual-null rep=", rr)
    Y0 <- make_residual_null_y(dat$Y, dat$X, cfg$seed + rr)
    for (mode in cfg$modes) {
      message("  mode=", mode)
      ii <- ii + 1L
      ans <- run_mode(Y0, dat$coords, dat$X, pathways, mode, rr, cfg)
      ans$dataset <- dataset
      rows[[ii]] <- ans
    }
  }
}
out <- do.call(rbind, rows)
out <- out[, c("dataset", "replicate", "mode", "pathway", "pathway_size", "p_H1",
               "p_H1_moment", "R2_cca", "PSVS_range", "mean_propSV", "elapsed_sec")]
write.csv(out, file.path(cfg$out_dir, "real_residual_null_pvalues.csv"), row.names = FALSE)
summary <- summarize_results(out)
write.csv(summary, file.path(cfg$out_dir, "real_residual_null_summary.csv"), row.names = FALSE)
md <- c(
  "# PASSAGE Real Residual Null Calibration",
  "",
  paste0("- Reps per dataset: ", cfg$n_reps),
  paste0("- Pathways per replicate: ", cfg$max_pathways),
  paste0("- Permutations: ", cfg$n_perm),
  "",
  paste(colnames(summary), collapse = " | "),
  paste(rep("---", ncol(summary)), collapse = " | "),
  apply(summary, 1L, function(x) paste(x, collapse = " | "))
)
writeLines(md, file.path(cfg$out_dir, "summary.md"))
message("Wrote outputs to ", cfg$out_dir)
print(summary, row.names = FALSE)
