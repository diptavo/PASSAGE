if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}

internal <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) get(name, mode = "function")
  else getFromNamespace(name, "PASSAGE")
}
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

set.seed(871L)
n <- 30L
coords <- cbind(x = runif(n), y = runif(n))
X <- cbind(intercept = 1, nuisance = rnorm(n))
Y <- matrix(rnorm(n * 10L), n, 10L)
Y[, 1:3] <- Y[, 1:3] + 2 * sin(coords[, 1L] * 4)
dimnames(Y) <- list(paste0("aggregation_spot_", seq_len(n)),
                    paste0("aggregation_gene_", seq_len(ncol(Y))))
rownames(coords) <- rownames(X) <- rownames(Y)
fit <- passage_fit_pathway(Y, colnames(Y), X = X, K = 3L)
cache <- passage_prepare_gp_kernel(coords, X = X, ranges = c(0.1, 0.3),
                                    kernel_method = "dense")
fractions <- c(0, 0.1, 0.4, 0.8)

# A statistic-only evaluation must neither consume nor initialize the RNG.
old_seed <- .Random.seed
details <- passage_pathway_gp_statistic(fit, cache, fractions, output = "details")
stopifnot(identical(.Random.seed, old_seed), is.finite(details$statistic),
          details$statistic >= 0, nrow(details$component_fits) == fit$K,
          identical(details$estimation$explained_variance, fit$explained_variance),
          all(c("range", "fraction", "statistic", "selected") %in%
                names(details$component_fits)))
rm(".Random.seed", envir = .GlobalEnv)
scalar <- passage_pathway_gp_statistic(fit, cache, fractions)
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
          identical(scalar, details$statistic))
assign(".Random.seed", old_seed, envir = .GlobalEnv)

contributions <- details$component_fits$statistic
stopifnot(abs(scalar - sum(contributions)) < 1e-12)
for (aggregation in c("sum", "max", "topk")) {
  expected <- switch(aggregation, sum = sum(contributions),
                     max = max(contributions),
                     topk = sum(sort(contributions, decreasing = TRUE)[1:2]))
  observed <- passage_pathway_gp_statistic(fit, cache, fractions,
                                           aggregation = aggregation,
                                           output = "details")
  calibrated <- passage_test_pathway_gp(fit, kernel_cache = cache,
      fractions = fractions, aggregation = aggregation, top_k = 2L,
      B = 23L, seed = 133L, batch_size = 5L)
  stopifnot(abs(observed$statistic - expected) < 1e-12,
            identical(calibrated$statistic, observed$statistic),
            identical(calibrated$component_fits, observed$component_fits),
            identical(calibrated$aggregation, aggregation),
            sum(observed$component_fits$selected) ==
              switch(aggregation, sum = 3L, max = 1L, topk = 2L))

  # Reconstruct each complete rotation independently, then reselect its own
  # largest components. Reusing the observed component selection would fail.
  score <- qr.resid(qr(cache$X), fit$scores)
  score <- sweep(score, 2L, sqrt(colSums(score^2)), "/")
  gram_chol <- chol(crossprod(score))
  set.seed(133L)
  manual_null <- replicate(23L, {
    rotation <- internal(".passage_gp_rotation")(ncol(cache$basis), cache$df, gram_chol)
    pp <- crossprod(cache$projection_transform, rotation)
    component <- internal(".passage_gp_profile")(pp^2, rep(1, 3L), cache,
                                                   fractions)$statistic
    switch(aggregation, sum = sum(component), max = max(component),
           topk = sum(sort(component, decreasing = TRUE)[1:2]))
  })
  stopifnot(max(abs(manual_null - calibrated$null)) < 1e-10,
            calibrated$p == (1 + sum(manual_null >= calibrated$statistic -
              1e-10 * max(1, abs(calibrated$statistic)))) / 24)
}

# Capping top_k at the available rank handles variance-based ranks uniformly.
sum_test <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = fractions,
                                    aggregation = "sum", B = 19L, seed = 47L)
top_all <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = fractions,
                                   aggregation = "topk", top_k = 5L, B = 19L, seed = 47L)
stopifnot(top_all$top_k == 5L, top_all$top_k_used == 3L,
          abs(top_all$statistic - sum_test$statistic) < 1e-12,
          max(abs(top_all$null - sum_test$null)) < 1e-12,
          identical(top_all$p, sum_test$p))
max_test <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = fractions,
                                    aggregation = "max", B = 19L, seed = 47L)
top_one <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = fractions,
                                   aggregation = "topk", top_k = 1L, B = 19L, seed = 47L)
stopifnot(identical(top_one$statistic, max_test$statistic),
          identical(top_one$null, max_test$null), identical(top_one$p, max_test$p))
one_fit <- passage_fit_pathway(Y, colnames(Y), X = X, K = 1L)
one_stat <- passage_pathway_gp_statistic(one_fit, cache, fractions,
                                         aggregation = "topk", output = "details")
stopifnot(one_stat$top_k_used == 1L, nrow(one_stat$component_fits) == 1L)

for (aggregation in c("sum", "max", "topk")) {
  tied <- passage_test_pathway_gp(fit, kernel_cache = cache, fractions = 0,
      aggregation = aggregation, B = 9L, seed = 61L)
  stopifnot(tied$statistic == 0, all(tied$null == 0), tied$p == 1)
}
expect_error(passage_pathway_gp_statistic(fit, cache, aggregation = "best"), "arg")
expect_error(passage_test_pathway_gp(fit, kernel_cache = cache,
                                    aggregation = "best", B = 3L), "arg")
for (k in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) {
  expect_error(passage_pathway_gp_statistic(fit, cache, top_k = k), "top_k")
}
expect_error(passage_pathway_gp_statistic(fit, list()), "invalid kernel_cache")
expect_error(passage_pathway_gp_statistic(fit, cache, fractions = c(0.2, 0.5)), "fractions")
bad_fit <- fit
bad_fit$scores_source <- "gp_predictions"
expect_error(passage_pathway_gp_statistic(bad_fit, cache), "unsmoothed")
bad_fit <- fit
bad_fit$X[, 2L] <- rev(bad_fit$X[, 2L])
expect_error(passage_pathway_gp_statistic(bad_fit, cache), "exactly match")

cat("GP statistic decoupling and sum/max/top-k null calibration contracts passed\n")
