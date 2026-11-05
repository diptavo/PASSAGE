if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}
internal <- function(nm) if (exists(nm, mode = "function")) get(nm) else getFromNamespace(nm, "PASSAGE")

set.seed(491)
coords <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
                               y = seq(0, 1, length.out = 6)))
X <- cbind(intercept = 1, technical = rnorm(nrow(coords)))
Y <- matrix(rnorm(nrow(coords) * 80), nrow(coords), 80)
Y <- sweep(Y, 2, seq(.2, 5, length.out = ncol(Y)), "*") + X[, 2] %o% seq(.2, 3, length.out = ncol(Y))
Y[, 1:10] <- Y[, 1:10] + 2 * coords[, 1]
colnames(Y) <- paste0("g", seq_len(ncol(Y)))
engine <- passage_fit_engine_pca(Y, coords, X = X, K = 3, m = 5,
                                 range_grid = c(.1, .3), scale = TRUE, verbose = FALSE)
pre <- passage_h_precompute(engine)
stopifnot(pre$qr$rank == qr(X)$rank)
residuals <- internal("passage_engine_residuals")(engine, Y, precomp = pre)
stopifnot(max(abs(residuals - engine$residuals)) < 1e-10)
bins <- passage_make_gene_bins(Y, n_mean_bins = 2, n_detect_bins = 2)
stats <- c("score_z", "score_Q", "score_robust_z", "score_sparse_topk_z", "score_lse_z",
           "score_factor_hsic", "score_coherence", "score_activity_hotspot", "score_pERSA_z",
           "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z",
           "mean_propSV_conditional", "cEPSV", "ePSV", "cwPVE_trace", "cwPVE_top",
           "mean_abs_spatial_cor")
fun <- passage_make_statistic(engine, Y, stats, gene_bins = bins)
P <- 1:8
stopifnot(identical(names(fun(P)), stats), all(is.finite(fun(P))),
          identical(fun(P), fun(colnames(Y)[P])))
stopifnot(inherits(try(fun(c(1, 1000)), silent = TRUE), "try-error"),
          inherits(try(passage_make_statistic(engine, Y, "missing"), silent = TRUE), "try-error"),
          inherits(try(passage_make_statistic(engine, Y, score_weigth = "equal"), silent = TRUE), "try-error"))

# Scores calculated from raw Y with the scaled engine equal scores in its
# explicitly transformed units. This catches mixing scaled loadings with raw Y.
reference_engine <- engine
reference_engine$gene_center <- rep(0, ncol(Y))
reference_engine$gene_scale <- rep(1, ncol(Y))
reference_engine$scale <- FALSE
reference_fun <- passage_make_statistic(reference_engine, residuals, stats,
                                        gene_bins = bins)
stopifnot(max(abs(fun(P) - reference_fun(P))) < 1e-7)
direct <- passage_pathway_score_stat(engine, Y, P, pre, output = "z")
stopifnot(abs(fun(P)[["score_z"]] - direct) < 1e-10)

# Batched scores have the same set order, orientation and values as individual
# evaluation. Non-equal factor weights intentionally do not use this fast path.
sets <- list(target = P, background1 = 12:19, background2 = 41:48)
one <- passage_make_statistic(engine, Y)
batched <- attr(one, "batch")(sets)
stopifnot(identical(dimnames(batched), list(names(sets), "score_z")),
          max(abs(batched[, 1] - vapply(sets, function(s) one(s)[[1]], numeric(1)))) < 1e-8,
          is.null(attr(passage_make_statistic(engine, Y, score_weight_scheme = "var"), "batch")))

# Prespecified options reach the underlying statistic and repeat for each set.
sparse <- passage_make_statistic(engine, Y, "score_sparse_topk_z", topk_grid = c(2L, 4L))
sparse_direct <- passage_sparse_topk_pathway_score_stat(engine, Y, P, pre,
                                                        topk_grid = c(2L, 4L))
stopifnot(abs(sparse(P)[[1]] - sparse_direct) < 1e-10)
legacy <- passage_conditional_competitive_test(engine, Y, P, gene_bins = bins,
                                                statistic = "score_z", B = 9, seed = 2)
stopifnot(abs(legacy$permutation$statistic - one(P)[[1]]) < 1e-10,
          is.finite(legacy$permutation$p), identical(legacy$statistic, "score_z"))
cat("Statistic interface and scaled-engine checks passed\n")
