if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}

expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
set.seed(10292026)
n <- 96L
X <- cbind(intercept = 1, technical = rnorm(n))
latent <- matrix(rnorm(n * 2L), n, 2L)
loadings <- rbind(c(1, 0.2), c(0.9, -0.1), c(1.1, 0.3), c(0.8, 0.1),
                  c(0.2, 1), c(-0.1, 0.9), c(0.3, 1.1), c(0.1, 0.8),
                  c(0.6, 0.6), c(0.7, -0.5))
Yp <- latent %*% t(loadings) + matrix(rnorm(n * 10L, sd = 0.35), n, 10L)
Yp <- Yp + X %*% matrix(rnorm(20L), 2L, 10L)
Y <- cbind(Yp, matrix(rnorm(n * 14L), n, 14L))
colnames(Y) <- paste0("g", seq_len(ncol(Y)))
ids <- paste0("spot", seq_len(n))
coords <- cbind(x = runif(n), y = runif(n))
rownames(Y) <- rownames(X) <- rownames(coords) <- ids
pathway <- colnames(Y)[seq_len(10L)]

# A pathway fit must not depend on genes outside the requested set.
fit <- passage_fit_pathway(Y, pathway, X, K = 3L, seed = 13L)
Y_other <- Y
Y_other[, 11:24] <- 1000 * Y_other[, 11:24] + coords[, 1L] * 1e5
fit_other <- passage_fit_pathway(Y_other, pathway, X, K = 3L, seed = 13L)
stopifnot(identical(fit$scores, fit_other$scores),
          identical(fit$loadings, fit_other$loadings),
          identical(fit$gene_names, pathway),
          identical(fit$scores_source, "unsmoothed_expression"))
tested <- passage_test_pathway_fit(fit, coords, B = 39L, seed = 19L)
tested_other <- passage_test_pathway_fit(fit_other, coords, B = 39L, seed = 19L)
stopifnot(identical(tested$statistic, tested_other$statistic),
          identical(tested$p, tested_other$p))

# The fitted scores/loadings reconstruct the covariate-adjusted truncated SVD.
res <- qr.resid(qr(X), Y[, pathway, drop = FALSE])
direct <- svd(res)
target <- direct$u[, 1:3] %*% diag(direct$d[1:3]) %*% t(direct$v[, 1:3])
stopifnot(max(abs(fit$scores %*% t(fit$loadings) - target)) < 1e-9,
          max(abs(crossprod(X, fit$scores))) < 1e-9,
          fit$residual_df == n - ncol(X))

# Explicit estimation followed by testing must reproduce the subspace API.
reference <- passage_subspace_test(Y[, pathway, drop = FALSE], coords, X = X,
                                  K = 3L, B = 39L, seed = 19L)
stopifnot(max(abs(tested$statistic - reference$statistic)) < 1e-8,
          identical(tested$p, reference$p),
          identical(tested$null, reference$null),
          abs(tested$explained_variance - fit$explained_variance) < 1e-12)

# Missing and constant genes are reported, while mismatched sample IDs are rejected.
Y_extra <- cbind(Y, constant_gene = 5)
with_missing <- passage_fit_pathway(Y_extra, c(pathway, "constant_gene", "absent_gene"),
                                    X, K = 3L)
stopifnot(identical(with_missing$missing_genes, "absent_gene"),
          identical(with_missing$dropped_genes, "constant_gene"),
          identical(with_missing$gene_names, pathway))
expect_error(passage_fit_pathway(Y, pathway, X[n:1, , drop = FALSE], K = 3L))
expect_error(passage_test_pathway_fit(fit, coords[n:1, , drop = FALSE], B = 9L))
kernel <- exp(-as.matrix(dist(coords)))
dimnames(kernel) <- list(rev(ids), ids)
expect_error(passage_test_pathway_fit(fit, spatial_kernel = kernel, B = 9L))

# Known singular values make the 95% rank choice unambiguous: two PCs suffice.
orth <- qr.Q(qr(cbind(X, matrix(rnorm(n * 4L), n, 4L))))[, 3:6, drop = FALSE]
rotation <- qr.Q(qr(matrix(rnorm(16L), 4L, 4L)))
known <- orth %*% diag(c(10, 5, 1, 0.1)) %*% t(rotation)
known <- known + X %*% matrix(rnorm(8L), 2L, 4L)
dimnames(known) <- list(ids, paste0("known", 1:4))
retained <- passage_fit_pathway(known, colnames(known), X, variance_threshold = 0.95)
stopifnot(retained$K == 2L,
          abs(retained$explained_variance - 125 / 126.01) < 1e-10)

# FA has its own fixed-rank rule and supplies finite unsmoothed regression scores.
expect_error(passage_fit_pathway(Y, pathway, X, method = "factor_analysis"))
fa <- passage_fit_pathway(Y, pathway, X, method = "factor_analysis", K = 2L, seed = 23L)
stopifnot(fa$K == 2L, is.na(fa$explained_variance),
          identical(fa$rank_rule, "prespecified rank"),
          all(is.finite(fa$scores)), all(is.finite(fa$loadings)),
          all(is.finite(fa$diagnostics$uniqueness)),
          all(is.finite(fa$diagnostics$noise_variance)))
loadings_table <- data.frame(gene = rownames(fa$loadings), fa$loadings)
stopifnot(nrow(loadings_table) == length(pathway), ncol(loadings_table) == fa$K + 1L)
fa_test <- passage_test_pathway_fit(fa, coords, B = 19L, seed = 29L)
stopifnot(all(is.finite(fa_test$p)), all(fa_test$p >= 1 / 20 & fa_test$p <= 1),
          identical(fa_test$estimation$method, "factor_analysis"))

# Descriptive GP outputs cannot accidentally enter the fixed-kernel test.
gp_object <- structure(list(scores = fit$scores), class = "passage_pathway_gp")
expect_error(passage_test_pathway_fit(gp_object, coords, B = 9L))
smoothed <- fit
smoothed$scores_source <- "gp_predictions"
expect_error(passage_test_pathway_fit(smoothed, coords, B = 9L))

# Fitting/testing preserve a caller's RNG state, including an initially absent seed.
set.seed(31L)
state <- .Random.seed
invisible(passage_fit_pathway(Y, pathway, X, method = "factor_analysis", K = 2L, seed = 37L))
stopifnot(identical(state, .Random.seed))
invisible(passage_test_pathway_fit(fit, coords, B = 9L, seed = 41L))
stopifnot(identical(state, .Random.seed))
rm(".Random.seed", envir = .GlobalEnv)
invisible(passage_fit_pathway(Y, pathway, X, K = 3L, seed = 43L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
invisible(passage_test_pathway_fit(fit, coords, B = 9L, seed = 47L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assign(".Random.seed", state, envir = .GlobalEnv)

cat("Pathway estimation, subspace equivalence, alignment, rank and RNG tests passed\n")
