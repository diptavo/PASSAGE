if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}

internal <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) get(name, mode = "function")
  else getFromNamespace(name, "PASSAGE")
}
expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
set.seed(10293001)
n <- 28L
ids <- paste0("inference_spot_", seq_len(n))
coords <- cbind(x = runif(n), y = runif(n))
X <- cbind(intercept = 1, nuisance = rnorm(n))
Y <- matrix(rnorm(n * 8L), n, 8L) + X %*% matrix(rnorm(16L), 2L, 8L)
dimnames(Y) <- list(ids, paste0("inference_gene_", seq_len(ncol(Y))))
rownames(coords) <- rownames(X) <- ids
fit <- passage_fit_pathway(Y, colnames(Y), X, K = 2L)
ranges <- c(0.15, 0.45)
fractions <- c(0, 0.25, 0.65, 0.9)
cache <- passage_prepare_gp_kernel(coords, X, ranges = ranges, kernel_method = "dense")
ans <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = fractions,
                                B = 39L, seed = 13L)
stopifnot(is.finite(ans$statistic), ans$statistic >= 0, length(ans$p) == 1L,
          is.finite(ans$p), ans$p >= 1 / 40, ans$p <= 1,
          length(ans$null) == 39L, all(is.finite(ans$null)),
          nrow(ans$component_fits) == 2L)

# Independently evaluate restricted/profile likelihood using dense inverses and
# determinants in an orthonormal residual basis. This checks the fast spectral
# calculation against the actual covariance model, not another implementation
# of its spectral formula.
qx <- qr(X)
df <- n - qx$rank
residual_basis <- qr.Q(qx, complete = TRUE)[, (qx$rank + 1L):n, drop = FALSE]
F <- crossprod(residual_basis, fit$scores)
ss <- colSums(F^2)
brute <- rep(0, ncol(F))
for (range in ranges) {
  z <- sqrt(3) * as.matrix(dist(coords)) / range
  R <- crossprod(residual_basis, ((1 + z) * exp(-z)) %*% residual_basis)
  R <- R * df / sum(diag(R))
  for (h in fractions) {
    C <- (1 - h) * diag(df) + h * R
    weighted_ss <- colSums(F * solve(C, F))
    improvement <- df * log(ss / weighted_ss) - as.numeric(determinant(C, logarithm = TRUE)$modulus)
    brute <- pmax(brute, improvement)
  }
}
stopifnot(max(abs(ans$component_fits$statistic - brute)) < 1e-7,
          abs(ans$statistic - sum(brute)) < 1e-7)
projection <- crossprod(cache$projection_transform, crossprod(cache$basis, fit$scores))
profile <- internal(".passage_gp_profile")(projection^2, colSums(fit$scores^2), cache, fractions)
stopifnot(max(abs(profile$statistic - brute)) < 1e-7)

# A grid containing only the null model has no improvement and all ties count
# toward the upper tail: plus-one Monte Carlo p must be exactly one.
null_only <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = 0,
                                      B = 19L, seed = 17L)
stopifnot(abs(null_only$statistic) < 1e-10, all(abs(null_only$null) < 1e-10),
          null_only$p == 1)
stopifnot(abs(ans$p - (1 + sum(ans$null >= ans$statistic)) / 40) < 1e-12)

# At full landmark rank, the PSD Nystrom construction must reproduce the dense
# covariance model (up to numerical tolerance) on the same observations.
full_lowrank <- passage_prepare_gp_kernel(coords, X, ranges = ranges,
                                           kernel_method = "lowrank", n_features = n)
full_lowrank_test <- passage_test_pathway_gp(fit, kernel_cache = full_lowrank,
                                             fractions = fractions, B = 19L, seed = 23L)
stopifnot(abs(full_lowrank_test$statistic - ans$statistic) < 1e-5)

# Profiling out scale must make the observed comparison invariant to arbitrary
# score units, component signs, and ordering. These changes can alter a specific
# Monte Carlo realization, so compare the observed statistic here.
rescaled <- fit
rescaled$scores <- sweep(fit$scores[, 2:1, drop = FALSE], 2L, c(-3.7, 0.24), "*")
rescaled_test <- passage_test_pathway_gp(rescaled, kernel_cache = cache,
                                         fractions = fractions, B = 9L, seed = 29L)
stopifnot(abs(rescaled_test$statistic - ans$statistic) < 1e-8)

# Check the joint rotation law with correlated FA-like scores. A whole-residual-
# space draw preserves the complete Gram, not merely each component's variance.
rotation <- internal(".passage_gp_rotation")
gram <- matrix(c(4, 1.8, 1.8, 2), 2L, 2L)
gram_chol <- chol(gram)
set.seed(31L)
whole <- rotation(q = 7L, df = 7L, gram_chol = gram_chol)
stopifnot(identical(dim(whole), c(7L, 2L)),
          max(abs(crossprod(whole) - gram)) < 1e-10)
near_whole <- rotation(q = 6L, df = 7L, gram_chol = gram_chol)
stopifnot(identical(dim(near_whole), c(6L, 2L)), all(is.finite(near_whole)),
          min(eigen(gram - crossprod(near_whole), symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
moments <- replicate(1200L, {
  draw <- rotation(q = 3L, df = 9L, gram_chol = gram_chol)
  c(draw[1L, ], as.vector(crossprod(draw)))
})
stopifnot(max(abs(rowMeans(moments)[1:2])) < 0.08,
          max(abs(matrix(rowMeans(moments)[3:6], 2L) - gram / 3)) < 0.12)

# PCA rank selection uses only the gene Gram. A row rotation fixing X changes
# scores by the same rotation, modulo arbitrary component signs, and preserves
# the rank and loading subspace. This supports conditioning rather than refitting
# covariance-only estimation in every null draw.
Q <- qr.Q(qr(matrix(rnorm(df * df), df, df)))
T <- residual_basis %*% Q %*% t(residual_basis) + X %*% solve(crossprod(X), t(X))
Y_rotated <- T %*% Y
dimnames(Y_rotated) <- dimnames(Y)
pca <- passage_fit_pathway(Y, colnames(Y), X, variance_threshold = 0.95)
pca_rotated <- passage_fit_pathway(Y_rotated, colnames(Y), X, variance_threshold = 0.95)
stopifnot(pca$K == pca_rotated$K,
          abs(pca$explained_variance - pca_rotated$explained_variance) < 1e-12,
          max(abs(tcrossprod(pca_rotated$scores) - T %*% tcrossprod(pca$scores) %*% t(T))) < 1e-8)
fa <- passage_fit_pathway(Y, colnames(Y), X, method = "factor_analysis", K = 1L)
fa_rotated <- passage_fit_pathway(Y_rotated, colnames(Y), X,
                                  method = "factor_analysis", K = 1L)
stopifnot(max(abs(tcrossprod(fa_rotated$scores) - T %*% tcrossprod(fa$scores) %*% t(T))) < 1e-6)
fa_test <- passage_test_pathway_gp(fa, kernel_cache = cache,
                                   fractions = fractions, B = 19L, seed = 35L)
stopifnot(is.finite(fa_test$p), fa_test$p >= 1 / 20, fa_test$p <= 1)

# Cached kernels are reusable only for the same observations and adjustment
# design. Invalid statistical controls and GP-smoothed input must fail clearly.
expect_error(passage_test_pathway_gp(fit, coords[n:1L, , drop = FALSE], B = 9L))
expect_error(passage_test_pathway_gp(fit, kernel_cache = cache, fractions = 1, B = 9L))
expect_error(passage_test_pathway_gp(fit, kernel_cache = cache, fractions = -0.1, B = 9L))
expect_error(passage_test_pathway_gp(fit, kernel_cache = cache, B = 0L))
expect_error(passage_prepare_gp_kernel(coords, X, ranges = 0))
bad_fit <- fit
bad_fit$scores_source <- "gp_predictions"
expect_error(passage_test_pathway_gp(bad_fit, kernel_cache = cache, B = 9L))
bad_fit <- fit
bad_fit$X[, 2L] <- rev(bad_fit$X[, 2L])
expect_error(passage_test_pathway_gp(bad_fit, kernel_cache = cache, B = 9L))

# Reproducibility must not consume a caller's RNG stream, including the case
# where that stream has not yet been initialized.
set.seed(37L)
saved <- .Random.seed
repeat_test <- passage_test_pathway_gp(fit, kernel_cache = cache,
                                        fractions = fractions, B = 39L, seed = 13L)
stopifnot(identical(saved, .Random.seed), identical(ans$null, repeat_test$null),
          identical(ans$p, repeat_test$p))
rm(".Random.seed", envir = .GlobalEnv)
invisible(passage_test_pathway_gp(fit, kernel_cache = cache, B = 9L, seed = 41L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assign(".Random.seed", saved, envir = .GlobalEnv)

cat("Integrated GP inference: covariance likelihood, null rotations, Nystrom equivalence, invariance and API contracts passed\n")
