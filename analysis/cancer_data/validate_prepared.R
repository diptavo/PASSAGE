#!/usr/bin/env Rscript
# Validate and smoke-test a complete exported sample with PASSAGE 0.6.1+.
# Rscript validate_prepared.R INPUT.rds REPORT.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply INPUT.rds REPORT.csv")
suppressPackageStartupMessages(library(PASSAGE))
suppressPackageStartupMessages(library(Matrix))
stopifnot(packageVersion("PASSAGE") >= "0.6.1")
b <- passage_validate_bundle(args[1L])
stopifnot(length(b$pathway_sets) == 50L, qr(b$X)$rank == ncol(b$X),
          identical(rownames(b$raw_counts), colnames(b$Y)),
          identical(colnames(b$raw_counts), rownames(b$Y)),
          identical(b$spot_data$barcode, rownames(b$Y)),
          all(is.finite(b$raw_counts@x)), all(b$raw_counts@x >= 0),
          all(b$raw_counts@x == round(b$raw_counts@x)))
coverage <- vapply(b$pathway_sets, function(g) sum(unique(g) %in% colnames(b$Y)), integer(1))
bins <- passage_make_gene_bins(b$Y, n_mean_bins = 3, n_detect_bins = 2, n_var_bins = 2)
supported <- vapply(b$pathway_sets, function(g) {
  passage_matching_diagnostics(g, bins, reference = "universe")$supported
}, logical(1))
stopifnot(all(coverage >= 15 & coverage <= 500), all(supported))

# Retain every spot and background gene. Low B checks execution, not significance.
started <- proc.time()[[3L]]
fit <- passage_analyze_bundle(
  b, pathways = b$pathway_sets["HALLMARK_HYPOXIA"], gene_bins = bins,
  B = 19, seed = 1, K = 3, m = 6, range_grid = c(0.1, 0.3),
  statistics = c("score_z", "score_activity_hotspot", "score_robust_z", "score_sparse_topk_z"),
  combine = "minp", reference = "universe", verbose = TRUE
)
stopifnot(nrow(fit$tests) == 5L, all(is.finite(fit$tests$p_competitive)),
          all(fit$tests$p_competitive > 0 & fit$tests$p_competitive <= 1))
xy <- sweep(b$coords, 2, apply(b$coords, 2, min), "-")
xy <- xy / max(xy)
pathway_fit <- passage_fit_pathway(b$Y, b$pathway_sets$HALLMARK_HYPOXIA,
                                  X = b$X, K = 3, seed = 1)
cache <- passage_prepare_gp_kernel(xy, X = b$X, ranges = c(0.1, 0.3), n_features = 32)
gp <- passage_test_pathway_gp(pathway_fit, kernel_cache = cache, B = 9, seed = 1)
stopifnot(is.finite(gp$p), gp$p > 0, gp$p <= 1)
report <- data.frame(sample = b$sample, n_spots = nrow(b$Y), n_genes = ncol(b$Y),
  n_covariates = ncol(b$X), hallmark_sets = length(coverage),
  hallmark_eligible = sum(coverage >= 15 & coverage <= 500),
  matching_supported = sum(supported), PASSAGE_version = as.character(packageVersion("PASSAGE")),
  md5 = unname(tools::md5sum(args[1L])), bytes = file.info(args[1L])$size,
  competitive_tests = nrow(fit$tests), competitive_B = 19, gp_B = 9,
  elapsed_seconds = proc.time()[[3L]] - started, status = "PASS")
dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
write.csv(report, args[2L], row.names = FALSE)
print(report)
