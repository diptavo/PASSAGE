# Run interactively after installing PASSAGE 0.6.0 or later.
# This small synthetic example is an interface demonstration.
library(PASSAGE)
set.seed(7)
coords <- as.matrix(expand.grid(x = seq(0, 1, length.out = 8),
                               y = seq(0, 1, length.out = 8)))
Y <- matrix(rnorm(64 * 160), 64, 160)
Y[, 1:12] <- Y[, 1:12] + 0.6 * sin(5 * coords[, 1])
dimnames(Y) <- list(paste0("spot", 1:64), paste0("gene", 1:160))
rownames(coords) <- rownames(Y)
pathway <- paste0("gene", 1:12)

bins <- passage_make_gene_bins(Y, n_mean_bins = 2,
                              n_detect_bins = 1, n_var_bins = 2)
support <- passage_matching_diagnostics(pathway, bins)
support$supported
support$strata

engine <- passage_fit_engine_pca(Y, coords, K = 3, m = 6,
                                range_grid = c(0.1, 0.3), verbose = FALSE)
statistic <- passage_make_statistic(
  engine, Y,
  statistics = c("score_z", "score_activity_hotspot", "score_sparse_topk_z"),
  gene_bins = bins
)
statistic(1:12)                 # numbers only; no p-values

result <- passage_competitive_test(
  pathway, statistic, gene_bins = bins,
  B = 99, seed = 1, combine = "cauchy"
)
result$components              # one competitive p-value per statistic
result$combinations            # calibrated combined competitive p-value
result$matching$strata         # available versus required genes in each bin
result$matching$balance        # observed versus background characteristics
