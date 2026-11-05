# PASSAGE 0.6.1: run this example in interactive R.
# Small B is for demonstrating the interface, not scientific inference.
library(PASSAGE)
set.seed(2026)
coords <- as.matrix(expand.grid(x = seq(0, 1, length.out = 8),
                               y = seq(0, 1, length.out = 8)))
X <- cbind(intercept = 1, technical = rnorm(64))
noise_sd <- exp(1.5 * coords[, 1])
Y <- matrix(rnorm(64 * 160), 64, 160) * noise_sd +
  X[, "technical"] %o% rnorm(160, sd = 2)
Y[, 1:12] <- Y[, 1:12] + sin(5 * coords[, 2])
dimnames(Y) <- list(paste0("spot", 1:64), paste0("gene", 1:160))
rownames(coords) <- rownames(X) <- rownames(Y)
pathway <- paste0("gene", 1:12)

# 1. Competitive: compare the pathway to full-universe matched gene sets.
engine <- passage_fit_engine_pca(Y, coords, X = X, K = 4,
  m = 8, range_grid = c(.1, .3), verbose = FALSE)
bins <- passage_make_gene_bins(Y, n_mean_bins = 2,
  n_detect_bins = 1, n_var_bins = 2)
statistic <- passage_make_statistic(engine, Y, gene_bins = bins,
  statistics = c("score_z", "score_activity_hotspot"))
competitive <- passage_competitive_test(pathway, statistic, bins,
  reference = "universe", combine = "minp", B = 19, seed = 1)
competitive$components
competitive$combinations

# 2. Self-contained: retain unequal residual magnitudes in bootstrap draws.
# Uses the measured mean design, with pathway estimation repeated every draw.
self_contained <- passage_test_pathway_gp_wild(Y, pathway, coords,
  X = X, K = 2, B = 19, seed = 1, n_features = 16)
self_contained$p
self_contained$component_fits

# 3. Optional negative-control design, if valid controls are known beforehand.
# These simulated genes contain technical variation but no added spatial signal.
# On real data, PASSAGE does not determine which genes are valid controls.
adjustment <- passage_prepare_nuisance(Y,
  control_genes = paste0("gene", 121:160), n_factors = 1)
adjustment$diagnostics
Y_without_controls <- Y[, adjustment$analysis_genes, drop = FALSE]
X_from_controls <- adjustment$X
# Rebuild bins, engine and statistic with Y_without_controls and X_from_controls.
# In self-contained inference, treating estimated factors as fixed requires
# additional validation; the bootstrap above used the known measured design.
