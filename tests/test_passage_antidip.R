set.seed(20260730)

root <- getwd()
if (dir.exists(file.path(root, "R"))) {
  for (f in sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 7),
  y = seq(0, 1, length.out = 7)
))
n <- nrow(coords)
g <- 42L
Y <- matrix(rnorm(n * g, sd = 0.8), nrow = n, ncol = g)
colnames(Y) <- paste0("g", seq_len(g))
X <- cbind(depth = as.numeric(scale(seq_len(n))))

d <- as.matrix(dist(coords))
z <- sqrt(3) * d / 0.22
Kc <- (1 + z) * exp(-z) + diag(1e-6, n)
u <- as.numeric(t(chol(Kc)) %*% rnorm(n))
for (jj in 1:8) {
  Y[, jj] <- Y[, jj] + 0.9 * u + rnorm(n, sd = 0.2)
}

pathways <- list(
  signal = paste0("g", 1:10),
  null = paste0("g", 21:32)
)
fit_args <- list(
  K = 2L,
  m = 5L,
  range_grid = c(0.10, 0.22, 0.40),
  kernel = "matern32",
  ordering = "none",
  verbose = FALSE
)

cf_engine <- passage_fit_engine_crossfit_spots(
  Y = Y,
  coords = coords,
  X = X,
  n_folds = 3L,
  seed = 1L,
  fit_args = fit_args,
  verbose = FALSE
)
stopifnot(inherits(cf_engine, "passage_engine"))
stopifnot(cf_engine$anti_dip == "spot_crossfit")
stopifnot(nrow(cf_engine$V) == n)
stopifnot(nrow(cf_engine$A) == g)

hold_engine <- passage_fit_engine_pathway_holdout(
  Y = Y,
  coords = coords,
  pathway = pathways$signal,
  X = X,
  fit_args = fit_args,
  verbose = FALSE
)
stopifnot(inherits(hold_engine, "passage_engine"))
stopifnot(hold_engine$anti_dip == "pathway_holdout")
stopifnot(!any(pathways$signal %in% hold_engine$background_genes))

pre_hold <- passage_h_precompute(hold_engine, X = X)
hold_score <- passage_score_test(
  hold_engine,
  Y,
  pathways$signal,
  pre_hold,
  gene_names = colnames(Y),
  calibration = "moment"
)
stopifnot(inherits(hold_score, "passage_score_result"))
stopifnot(is.finite(hold_score$p_omnibus))

refit_score <- passage_score_test_refit_null(
  Y = Y,
  coords = coords,
  pathway = pathways$signal,
  X = X,
  fit_args = fit_args,
  n_perm = 5L,
  seed = 2L,
  return_engine = TRUE,
  verbose = FALSE
)
stopifnot(inherits(refit_score, "passage_score_result"))
stopifnot(refit_score$anti_dip == "refit_null")
stopifnot(refit_score$calibration$method == "refit_null_residual_permutation")
stopifnot(length(refit_score$calibration$null_stat) == 5L)
stopifnot(inherits(refit_score$engine, "passage_engine"))

run_cf <- passage_run(
  Y = Y,
  coords = coords,
  X = X,
  pathways = pathways,
  K = 2L,
  m = 5L,
  range_grid = c(0.10, 0.22, 0.40),
  hypotheses = "H1",
  calibration = "moment",
  anti_dip = "spot_crossfit",
  anti_dip_args = list(n_folds = 3L, seed = 3L, fit_args = fit_args),
  verbose = FALSE
)
stopifnot(inherits(run_cf, "passage_run"))
stopifnot(run_cf$anti_dip == "spot_crossfit")
stopifnot(all(run_cf$summary$status == "tested"))

run_hold <- passage_run(
  Y = Y,
  coords = coords,
  X = X,
  pathways = pathways,
  K = 2L,
  m = 5L,
  range_grid = c(0.10, 0.22, 0.40),
  hypotheses = "H1",
  calibration = "moment",
  anti_dip = "pathway_holdout",
  anti_dip_args = list(fit_args = fit_args),
  verbose = FALSE
)
stopifnot(inherits(run_hold, "passage_run"))
stopifnot(run_hold$anti_dip == "pathway_holdout")
stopifnot(all(run_hold$summary$status == "tested"))

run_refit <- passage_run(
  Y = Y,
  coords = coords,
  X = X,
  pathways = pathways[1],
  K = 2L,
  m = 5L,
  range_grid = c(0.10, 0.22, 0.40),
  hypotheses = "H1",
  calibration = "moment",
  n_perm = 5L,
  anti_dip = "refit_null",
  anti_dip_args = list(fit_args = fit_args, n_perm = 5L),
  verbose = FALSE
)
stopifnot(inherits(run_refit, "passage_run"))
stopifnot(run_refit$anti_dip == "refit_null")
stopifnot(is.finite(run_refit$summary$p_H1))

cat("PASSAGE anti-double-dipping tests passed\n")
