if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else library(PASSAGE)
expect_error <- function(x) stopifnot(inherits(try(force(x), silent = TRUE), "try-error"))

# Exhaust every possible label set. Full-universe references are independent of
# the observed labels; finite-population tail ranks must be super-uniform.
genes <- paste0("g", 1:8)
sets <- combn(1:8, 3, simplify = FALSE)
score <- function(idx) c(total = sum(c(1, 2, 2, 4, 7, 9, 12, 15)[idx]), peak = max(idx))
p <- vapply(sets, function(P) {
  a <- passage_competitive_test(P, score, gene_names = genes,
    matched_sets = sets, reference = "universe", combine = c("cauchy", "minp", "fisher"))
  c(a$components$p_competitive, a$combinations$p_competitive)
}, numeric(5))
for (alpha in c(.05, .1, .25, .5)) stopifnot(all(rowMeans(p <= alpha) <= alpha + 1e-12))
expect_error(passage_competitive_test(1:3, score, gene_names = genes,
  matched_sets = sets, reference = "complement"))
one <- data.frame(gene = genes, bin = seq_along(genes))
a <- passage_competitive_test(1:3, score, gene_bins = one, reference = "universe", B = 9)
stopifnot(all(a$components$p_competitive == 1), !a$matching$reference_has_variation,
  all(vapply(a$matched_sets, function(s) setequal(s, 1:3), logical(1))))
expect_error(passage_competitive_test(1:3, score, gene_bins = one, B = 9))
set.seed(60); old <- .Random.seed
a <- passage_competitive_test(1:3, score, gene_names = genes, B = 49, reference = "universe")
stopifnot(identical(old, .Random.seed), a$matching$mean_target_overlap > 0,
  all(vapply(a$matched_sets, function(s) !anyDuplicated(s), logical(1))))

# Covariates and control factors: require alignment, keep the measured design,
# exclude controls from the reported analysis universe, and report leverage.
n <- 36L
Y <- matrix(rnorm(n * 30), n, 30, dimnames = list(paste0("s", 1:n), paste0("g", 1:30)))
X <- cbind(intercept = 1, measured = rnorm(n)); rownames(X) <- rownames(Y)
Y <- Y + X[, 2] %o% seq_len(30) / 10
adjust <- passage_prepare_nuisance(Y, X, colnames(Y)[21:30], n_factors = 2)
stopifnot(ncol(adjust$X) == 4, qr(adjust$X)$rank == 4,
  identical(adjust$X[, 1:2], X), length(adjust$analysis_genes) == 20,
  max(abs(crossprod(X, adjust$X[, 3:4]))) < 1e-9,
  abs(sum(adjust$diagnostics$leverage) - 4) < 1e-10)
Ychanged <- Y; Ychanged[, 1:20] <- Ychanged[, 1:20] + 100
adjust2 <- passage_prepare_nuisance(Ychanged, X, colnames(Y)[21:30], 2)
stopifnot(identical(adjust$X, adjust2$X))
expect_error(passage_prepare_nuisance(Y, X, n_factors = 2))
expect_error(passage_prepare_nuisance(Y, X, c("g1", "absent"), 1))

# The robust test changes calibration only: the observed GP statistic must
# equal the original. Null draws must refit exactly the same PCA and design.
coords <- cbind(x = rep(seq(0, 1, length.out = 6), 6), y = rep(seq(0, 1, length.out = 6), each = 6))
rownames(coords) <- rownames(Y)
fit <- passage_fit_pathway(Y, 1:12, X = X, K = 2)
cache <- passage_prepare_gp_kernel(coords, X = X, ranges = c(.1, .3), kernel_method = "dense")
old <- .Random.seed
a <- passage_test_pathway_gp_wild(Y, 1:12, X = X, K = 2, kernel_cache = cache,
  B = 7, seed = 42, fractions = c(0, .2, .7))
stopifnot(identical(.Random.seed, old), all(a$null_ranks == 2),
  identical(a$statistic, passage_pathway_gp_statistic(fit, cache, fractions = c(0, .2, .7))))
q <- qr(X); r <- qr.resid(q, Y[, 1:12]); f <- Y[, 1:12] - r
h <- rowSums(qr.Q(q)^2); r <- r / sqrt(1 - h)
set.seed(42)
manual <- replicate(7, {
  w <- 2 * sample.int(2, n, replace = TRUE) - 3
  bf <- passage_fit_pathway(f + r * w, colnames(r), X = X, K = 2, seed = 42)
  passage_pathway_gp_statistic(bf, cache, fractions = c(0, .2, .7))
})
stopifnot(max(abs(manual - a$null)) < 1e-10)
tie <- passage_test_pathway_gp_wild(Y, 1:12, X = X, K = 2, kernel_cache = cache,
  B = 3, fractions = 0, keep_null = FALSE)
stopifnot(tie$p == 1, is.null(tie$null), is.null(tie$null_ranks))
fa <- passage_test_pathway_gp_wild(Y, 1:12, X = X, K = 1, method = "factor_analysis",
  kernel_cache = cache, B = 3)
stopifnot(is.finite(fa$p), all(fa$null_ranks == 1))
expect_error(passage_test_pathway_gp_wild(Y, 1:12, X = X, K = 2, kernel_cache = cache, ranges = .7, B = 3))
cat("Full-universe randomization, negative-control adjustment and wild-refit contracts passed\n")
