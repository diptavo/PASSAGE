#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: passage_pathway_cli.R BUNDLE.rds PATHWAY_NAME OUTPUT_DIR",
  "[B=9999] [SEED=1] [METHOD=pca] [K=auto] [RANGE=0.25] [DESCRIPTIVE_GP=false]",
  "[TEST=subspace] [GP_KERNEL=lowrank] [FEATURES=64]",
  "[RANGES=0.05,0.1,0.2,0.4] [FRACTIONS=0,0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.97]")
if (length(args) < 3L || length(args) > 14L) stop(usage)
suppressPackageStartupMessages(library(PASSAGE))
integer_arg <- function(position, default, label, minimum = 1L) {
  if (length(args) < position) return(default)
  x <- suppressWarnings(as.numeric(args[position]))
  if (length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < minimum || x > .Machine$integer.max)
    stop(label, " must be an integer >= ", minimum)
  as.integer(x)
}
grid_arg <- function(position, default, label) {
  if (length(args) < position) return(default)
  x <- suppressWarnings(as.numeric(strsplit(args[position], ",", fixed = TRUE)[[1L]]))
  if (!length(x) || any(!is.finite(x))) stop(label, " must be comma-separated finite numbers")
  x
}
B <- integer_arg(4L, 9999L, "B")
seed <- integer_arg(5L, 1L, "SEED", minimum = 0L)
method <- if (length(args) >= 6L) args[6L] else "pca"
if (!method %in% c("pca", "factor_analysis")) stop("METHOD must be pca or factor_analysis")
K <- if (length(args) >= 7L && args[7L] != "auto") integer_arg(7L, NULL, "K") else NULL
rho <- if (length(args) >= 8L) suppressWarnings(as.numeric(args[8L])) else 0.25
if (!is.finite(rho) || rho <= 0) stop("RANGE must be positive and finite")
gp_flag <- if (length(args) >= 9L) tolower(args[9L]) else "false"
if (!gp_flag %in% c("true", "false")) stop("DESCRIPTIVE_GP must be true or false")
test_choice <- if (length(args) >= 10L) args[10L] else "subspace"
if (!test_choice %in% c("subspace", "gp")) stop("TEST must be subspace or gp")
gp_kernel <- if (length(args) >= 11L) args[11L] else "lowrank"
if (!gp_kernel %in% c("lowrank", "dense")) stop("GP_KERNEL must be lowrank or dense")
n_features <- integer_arg(12L, 64L, "FEATURES")
ranges <- grid_arg(13L, c(0.05, 0.1, 0.2, 0.4), "RANGES")
fractions <- grid_arg(14L, c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97), "FRACTIONS")
if (any(ranges <= 0)) stop("RANGES must be positive")
if (any(fractions < 0 | fractions >= 1) || !0 %in% fractions)
  stop("FRACTIONS must include 0 and lie in [0, 1)")
bundle <- passage_validate_bundle(args[1L])
pathway <- args[2L]
if (is.null(bundle$pathway_sets[[pathway]])) stop("pathway not in bundle$pathway_sets")
out <- args[3L]
# Never leave an old completion marker beside a failed rerun.
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("OUTPUT_DIR must be new or empty; use a separate directory for each run")
xy <- sweep(bundle$coords, 2L, apply(bundle$coords, 2L, min), "-")
span <- max(xy)
if (!is.finite(span) || span <= 0) stop("coordinates have zero extent")
xy <- xy / span
fit <- passage_fit_pathway(bundle$Y, bundle$pathway_sets[[pathway]],
  X = bundle$X, method = method, K = K, seed = seed)
if (test_choice == "subspace") {
  test <- passage_test_pathway_fit(fit, xy, kernel_range = rho, B = B, seed = seed)
} else {
  test <- passage_test_pathway_gp(fit, coords = xy, ranges = ranges,
    fractions = fractions, B = B, seed = seed, kernel_method = gp_kernel,
    n_features = n_features)
}
gp <- if (gp_flag == "true") passage_fit_pathway_gp(fit, xy, seed = seed) else NULL
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(out, "pathway_fit.rds"))
saveRDS(test, file.path(out, "pathway_test.rds"))
write.csv(data.frame(spot = rownames(fit$scores), fit$scores),
  file.path(out, "component_scores.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(fit$loadings), fit$loadings),
  file.path(out, "component_loadings.csv"), row.names = FALSE)
summary <- data.frame(pathway = pathway, method = method, K = fit$K,
  n_spots = nrow(fit$scores), n_genes = length(fit$gene_names), B = B,
  trace = NA_real_, largest_root = NA_real_, p_trace = NA_real_, p_root = NA_real_,
  p_omnibus = NA_real_, statistic_gp = NA_real_, p_gp = NA_real_,
  testing_method = test_choice,
  inference = "self-contained Gaussian reference; no across-pathway multiplicity correction")
if (test_choice == "subspace") {
  summary$trace <- unname(test$statistic["trace"])
  summary$largest_root <- unname(test$statistic["largest_root"])
  summary$p_trace <- unname(test$p["trace"])
  summary$p_root <- unname(test$p["largest_root"])
  summary$p_omnibus <- test$p_omnibus_bonferroni
} else {
  summary$statistic_gp <- test$statistic
  summary$p_gp <- test$p
  write.csv(test$component_fits, file.path(out, "gp_component_fits.csv"), row.names = FALSE)
}
write.csv(summary, file.path(out, "pathway_summary.csv"), row.names = FALSE)
if (!is.null(gp)) {
  saveRDS(gp, file.path(out, "descriptive_gp.rds"))
  write.csv(gp$summary, file.path(out, "descriptive_gp_summary.csv"), row.names = FALSE)
}
saveRDS(list(sample = bundle$sample, pathway = pathway,
  version = as.character(packageVersion("PASSAGE")), B = B, seed = seed,
  estimation_method = method, rank_requested = K, rank_retained = fit$K,
  testing_method = test_choice, descriptive_gp = gp_flag == "true",
  coordinate_offset = apply(bundle$coords, 2L, min), coordinate_scale = span,
  kernel_range = if (test_choice == "subspace") rho else NULL,
  gp_configuration = if (test_choice == "gp") test$configuration else NULL,
  source = normalizePath(args[1L]), source_md5 = tools::md5sum(args[1L])),
  file.path(out, "configuration.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
writeLines("PASS: execution completed; scientific calibration remains required", file.path(out, "analysis.complete"))
