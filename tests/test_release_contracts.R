if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
  internal <- function(x) get(x, mode = "function")
} else {
  library(PASSAGE)
  internal <- function(x) getFromNamespace(x, "PASSAGE")
}
set.seed(62022)
Y <- scale(matrix(rnorm(80 * 40), 80, 40), scale = FALSE)
s <- svd(Y)
for (method in c("pca", "varimax")) {
  f <- internal("passage_low_rank_basis")(Y, s, 5L, s$d[1:5]^2 / 79, method)
  target <- s$u[, 1:5] %*% diag(s$d[1:5]) %*% t(s$v[, 1:5])
  stopifnot(max(abs(f$V %*% t(f$A) - target)) < 1e-10)
  stopifnot(max(abs(crossprod(f$V) / 79 - diag(5))) < 1e-10)
}
# Singleton candidate pools must never invoke sample(positive_integer).
bins <- data.frame(gene = letters[1:5], bin = c("a", "b", "c", "d", "a"))
draws <- internal("passage_sample_matched_indices")(1L, bins, B = 40, seed = 1)
stopifnot(all(vapply(draws, identical, logical(1), 5L)))
sets <- passage_sample_unique_sets(c(1, 2, 6), rep(c("a", "b"), each = 5), 30, 12)
stopifnot(all(vapply(sets, function(x) length(x) == 3 && !anyDuplicated(x) &&
  !any(x %in% c(1, 2, 6)) && sum(x <= 5) == 2, logical(1))))
stopifnot(inherits(try(passage_sample_unique_sets(1:4, rep("a", 5), 10), silent = TRUE), "try-error"))
coords <- cbind(runif(80), runif(80))
rownames(Y) <- rownames(coords) <- paste0("s", 1:80)
colnames(Y) <- paste0("g", 1:40)
obj <- list(Y = Y, coords = coords, sample = "synthetic")
stopifnot(nrow(passage_validate_bundle(obj)$X) == 80)
bad <- obj; bad$coords <- bad$coords[80:1, ]
stopifnot(inherits(try(passage_validate_bundle(bad), silent = TRUE), "try-error"))
# Trace/root statistics equal the explicit generalized quadratic forms.
ans <- passage_subspace_test(Y, coords, K = 5, B = 39, seed = 1)
z <- sqrt(3) * as.matrix(dist(coords)) / 0.25
R <- (1 + z) * exp(-z)
F <- s$u[,1:5] %*% diag(s$d[1:5]); A <- crossprod(F, R %*% F); C <- crossprod(F)
direct <- solve(C, A)
stopifnot(abs(ans$statistic["trace"] - 79 * sum(diag(direct))) < 1e-8)
stopifnot(abs(ans$statistic["largest_root"] - max(Re(eigen(direct)$values))) < 1e-8)
scaled <- passage_subspace_test(Y * 17, coords, K = 5, B = 39, seed = 1)
stopifnot(max(abs(ans$statistic - scaled$statistic)) < 1e-8)
stopifnot(identical(ans$p, scaled$p), all(ans$p >= 1/40 & ans$p <= 1))
eng <- passage_fit_engine_pca(Y, coords, K = 3, m = 5,
  range_grid = c(.1, .3), verbose = FALSE)
pre <- passage_h_precompute(eng)
sets <- list(1:10, 11:25, c(2,4,8,12,16))
rr <- internal("passage_residualize_with_qr")(Y, pre$qr)
fast <- internal("passage_score_set_batch")(eng, rr, sets, pre, chunk_size = 2L)
slow <- vapply(sets, function(p) passage_pathway_score_stat(eng, Y, p, pre), numeric(1))
stopifnot(max(abs(fast - slow)) < 1e-8)
cat("Reconstruction, unique matching, alignment and subspace tests passed\n")
