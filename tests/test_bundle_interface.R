if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}
close <- function(a, b, tol = 1e-8) {
  stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol,
                            check.attributes = FALSE)))
}
expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))

set.seed(69013L)
coords <- as.matrix(expand.grid(x = seq(2, 5, length.out = 6L),
                                y = seq(-1, 1, length.out = 6L)))
n <- nrow(coords)
X <- cbind(intercept = 1, technical = rnorm(n))
Y <- matrix(rnorm(n * 90L), n, 90L)
Y <- sweep(Y, 2L, seq(0.3, 3, length.out = ncol(Y)), "*") +
  X[, 2L] %o% rnorm(ncol(Y))
Y[, 1:8] <- Y[, 1:8] + sin(coords[, 1L])
dimnames(Y) <- list(paste0("bundle_spot_", seq_len(n)),
                    paste0("bundle_gene_", seq_len(ncol(Y))))
rownames(coords) <- rownames(X) <- rownames(Y)
pathways <- list(first = colnames(Y)[1:8], second = colnames(Y)[15:24])
bundle <- list(Y = Y, coords = coords, X = X, pathway_sets = pathways,
               sample = "interface_fixture")
# Fixed strata keep the fixture focused on interface behavior; matching itself
# has separate margin/support tests. Include real-valued balance diagnostics.
bins <- data.frame(gene = colnames(Y), bin = "all", mean_expr = colMeans(Y),
                    detection = colMeans(Y > 0), variance = apply(Y, 2L, var))
B <- 11L
seed <- 198L
saved_seed <- .Random.seed
default <- passage_analyze_bundle(bundle, B = B, seed = seed, K = 3L, m = 5L,
    range_grid = c(.1, .3), gene_bins = bins, min_size = 5L, verbose = FALSE,
    keep_matched_sets = TRUE)
stopifnot(identical(.Random.seed, saved_seed), nrow(default$summary) == 2L,
          nrow(default$tests) == 2L, all(default$summary$primary_test == "score_z"),
          default$configuration$estimation == "global", length(default$estimation) == 1L,
          default$estimation[[1L]]$design_rank == 2L,
          is.finite(default$estimation[[1L]]$rank_info$cumulative_variance))
compact <- passage_analyze_bundle(bundle, B = B, seed = seed, K = 3L, m = 5L,
    range_grid = c(.1, .3), gene_bins = bins, min_size = 5L, verbose = FALSE)
close(compact$summary, default$summary)
universe <- passage_analyze_bundle(bundle, B = B, seed = seed, K = 3L, m = 5L,
    range_grid = c(.1, .3), gene_bins = bins, min_size = 5L, verbose = FALSE,
    reference = "universe", keep_matched_sets = TRUE)
stopifnot(universe$configuration$reference == "universe",
          all(vapply(universe$details, function(d)
            d$competitive$matching$reference == "universe", logical(1))))
for (j in seq_along(pathways)) {
  expected_sets <- passage_sample_unique_sets(match(pathways[[j]], colnames(Y)),
    bins$bin, B, seed + j, reference = "universe")
  stopifnot(identical(universe$details[[j]]$competitive$matched_sets, expected_sets))
}
for (j in seq_along(pathways)) {
  stopifnot(is.null(compact$details[[j]]$competitive$matched_sets))
  close(compact$details[[j]]$competitive$scores, default$details[[j]]$competitive$scores)
}

# Independently perform the previous single-statistic workflow with the same
# seeded global fit and matching draws. Coordinate scaling must preserve aspect
# ratio and agree with the bundle contract.
scaled_coords <- sweep(coords, 2L, apply(coords, 2L, min), "-")
scaled_coords <- scaled_coords / max(scaled_coords)
set.seed(seed)
engine <- passage_fit_engine_pca(Y, scaled_coords, X = X, K = 3L, m = 5L,
                                  range_grid = c(.1, .3), verbose = FALSE)
pre <- passage_h_precompute(engine, X = X)
for (j in seq_along(pathways)) {
  P <- match(pathways[[j]], colnames(Y))
  sets <- passage_sample_unique_sets(P, bins$bin, B, seed + j)
  expected <- vapply(c(list(P), sets), function(idx)
    passage_pathway_score_stat(engine, Y, idx, pre, output = "z"), numeric(1))
  got <- default$details[[j]]$competitive
  close(got$scores[, 1L], expected)
  close(default$summary$score_z[j], expected[1L])
  close(default$summary$p_competitive[j], (1 + sum(expected[-1L] >= expected[1L])) / (B + 1))
  stopifnot(identical(got$matched_sets, sets), nrow(got$matching$balance) == 3L,
            got$matching$n_draws == B, is.finite(got$components$mc_lower),
            is.finite(got$components$mc_upper))
}

# Complementary statistics and every requested combiner share the same sets;
# primary output follows the explicitly requested combination ordering.
four <- c("score_z", "score_activity_hotspot", "score_robust_z", "score_sparse_topk_z")
combined <- passage_analyze_bundle(bundle, B = B, seed = seed, K = 3L, m = 5L,
    range_grid = c(.1, .3), gene_bins = bins, min_size = 5L, verbose = FALSE,
    statistics = four, combine = c("minp", "cauchy", "fisher"), keep_matched_sets = TRUE)
stopifnot(nrow(combined$tests) == length(pathways) * 7L,
          all(combined$summary$primary_test == "combined_minp"),
          all(is.na(combined$summary$mc_lower)), all(is.na(combined$summary$mc_upper)))
for (j in seq_along(pathways)) {
  got <- combined$details[[j]]$competitive
  stopifnot(identical(colnames(got$scores), four), all(is.finite(got$scores)),
            identical(got$matched_sets, default$details[[j]]$competitive$matched_sets),
            identical(got$combinations$method, c("minp", "cauchy", "fisher")))
  close(got$scores[, "score_z"], default$details[[j]]$competitive$scores[, "score_z"])
  close(combined$summary$p_competitive[j], got$combinations$p_competitive[1L])
}
for (label in unique(combined$tests$test)) {
  take <- combined$tests$test == label
  close(combined$tests$q_BH[take], p.adjust(combined$tests$p_competitive[take], "BH"))
  close(combined$tests$q_BY[take], p.adjust(combined$tests$p_competitive[take], "BY"))
}

# A prespecified rank grid produces separate evidence components, not an
# uncorrected best-rank p-value. Technical adjustment and scaling remain active.
ranked <- passage_analyze_bundle(bundle, B = 5L, seed = seed, K = 3L, m = 5L,
    range_grid = c(.1, .3), gene_bins = bins, min_size = 5L, verbose = FALSE,
    statistics = c("score_z", "score_activity_hotspot"), combine = "minp",
    rank_grid = c(2L, 3L), scale = TRUE)
stopifnot(identical(names(ranked$estimation), c("K2", "K3")),
          all(vapply(ranked$estimation, function(z) isTRUE(z$scaled_genes), logical(1))),
          all(vapply(ranked$estimation, function(z) z$design_rank == 2L, logical(1))),
          all(is.finite(ranked$tests$p_competitive)),
          all(is.na(ranked$summary$score_z)))
expected_names <- c("score_z_K2", "score_activity_hotspot_K2",
                     "score_z_K3", "score_activity_hotspot_K3")
stopifnot(identical(colnames(ranked$details[[1L]]$competitive$scores), expected_names))
expect_error(passage_analyze_bundle(bundle, B = 3L, gene_bins = bins, min_size = 5L,
                                    verbose = FALSE, rank_grid = c(2L, 2L)))

# Pathway-local GP fits are repeated once per target/background set and the
# profile contributions are reused across all three aggregation choices.
local <- passage_analyze_bundle(bundle, B = 5L, seed = seed, gene_bins = bins,
    min_size = 5L, verbose = FALSE, estimation = "pathway", local_K = 3L,
    gp_aggregations = c("sum", "max", "topk"), gp_top_k = 2L,
    gp_kernel_method = "dense", gp_ranges = c(.1, .3), combine = "cauchy",
    keep_matched_sets = TRUE)
local_names <- c("local_gp_sum", "local_gp_max", "local_gp_topk2")
stopifnot(local$configuration$estimation == "pathway",
          all(local$summary$primary_test == "combined_cauchy"),
          all(is.na(local$summary$score_z)), all(is.na(local$summary$ePSV)),
          all(local$summary$candidate_drivers == ""))
cache <- passage_prepare_gp_kernel(scaled_coords, X = X, ranges = c(.1, .3),
                                    kernel_method = "dense")
for (j in seq_along(pathways)) {
  got <- local$details[[j]]$competitive
  stopifnot(identical(colnames(got$scores), local_names),
            length(got$estimation) == 6L,
            all(vapply(got$estimation, function(z) z$K == 3L, logical(1))),
            nrow(got$estimation[[1L]]$component_fits) == 3L,
            all(vapply(got$estimation[-1L], function(z) is.null(z$component_fits), logical(1))))
  sets <- c(list(match(pathways[[j]], colnames(Y))), got$matched_sets)
  expected <- t(vapply(sets, function(idx) {
    fit <- passage_fit_pathway(Y, idx, X = X, K = 3L)
    detail <- passage_pathway_gp_statistic(fit, cache, output = "details")
    component <- detail$component_fits$statistic
    c(sum(component), max(component), sum(sort(component, decreasing = TRUE)[1:2]))
  }, numeric(3L)))
  close(got$scores, expected)
}

# Variance-based local ranks are recorded for each refitted background, with no
# rank borrowing from the observed pathway. The diagnostics reset per test.
factory <- passage_make_local_gp_statistic(Y, scaled_coords, X = X,
    kernel_cache = cache, variance_threshold = .80,
    aggregations = c("sum", "max", "topk"), top_k = 2L)
local_test <- passage_competitive_test(1:8, factory, gene_bins = bins, B = 3L, seed = seed)
stopifnot(length(local_test$estimation) == 4L,
          all(vapply(local_test$estimation, function(z) z$explained_variance >= .80, logical(1))))
local_repeat <- passage_competitive_test(15:24, factory, gene_bins = bins, B = 2L, seed = seed)
stopifnot(length(local_repeat$estimation) == 3L)

# Labeled covariates must never be silently relabeled into a different order.
misaligned_X <- X[n:1L, , drop = FALSE]
expect_error(passage_make_local_gp_statistic(Y, scaled_coords, X = misaligned_X,
                                             K = 3L))
expect_error(passage_make_local_gp_statistic(Y, scaled_coords[n:1L, , drop = FALSE],
                                             X = X, K = 3L))
# A cached kernel is reusable only for the same geometry, adjustment design,
# and any explicitly supplied kernel controls. Matching IDs alone is not enough.
other_coords <- scaled_coords
other_coords[1L, 1L] <- other_coords[1L, 1L] + 0.1
expect_error(passage_make_local_gp_statistic(Y, other_coords, X = X, kernel_cache = cache))
other_X <- X
other_X[, 2L] <- rev(other_X[, 2L])
expect_error(passage_make_local_gp_statistic(Y, scaled_coords, X = other_X,
                                             kernel_cache = cache))
expect_error(passage_make_local_gp_statistic(Y, scaled_coords, X = X,
                                             kernel_cache = cache, ranges = c(.1, .4)))
expect_error(passage_make_local_gp_statistic(Y, scaled_coords, X = X,
                                             kernel_cache = cache, kernel_method = "lowrank"))
expect_error(passage_make_local_gp_statistic(Y, scaled_coords, X = X,
                                             kernel_cache = cache, n_features = 65L))
compatible <- passage_make_local_gp_statistic(Y, scaled_coords, X = X, K = 3L,
    kernel_cache = cache, ranges = c(.3, .1), kernel_method = "dense", n_features = 64L)
stopifnot(all(is.finite(compatible(1:8))))

# Public matching and bundle boundaries must reject invalid numeric indices
# before a resolver can truncate fractional values or drop out-of-range genes.
for (bad in list(c(1, 2.5), c(0, 1), c(1, ncol(Y) + 1L), c(1, NA_real_), c(1, Inf))) {
  expect_error(passage_matching_diagnostics(bad, bins))
  expect_error(passage_competitive_test(bad, function(idx) sum(idx),
                                        gene_bins = bins, B = 3L))
}
expect_error(passage_analyze_bundle(bundle, pathways = list(bad = c(1:8, 9.5)),
    B = 3L, gene_bins = bins, min_size = 5L, verbose = FALSE))
expect_error(passage_analyze_bundle(bundle, pathways = list(bad = c(1:8, ncol(Y) + 1L)),
    B = 3L, gene_bins = bins, min_size = 5L, verbose = FALSE))

# A user callback does not require any PASSAGE fitting engine. It is evaluated
# on precisely the observed and matched gene-index sets.
gene_score <- colMeans(Y^2)
custom <- function(idx) c(total = sum(gene_score[idx]), peak = max(gene_score[idx]))
custom_result <- passage_analyze_bundle(bundle, B = 5L, seed = seed, gene_bins = bins,
    min_size = 5L, verbose = FALSE, statistics = custom, combine = "fisher",
    keep_matched_sets = TRUE)
stopifnot(custom_result$configuration$estimation == "custom",
          length(custom_result$estimation) == 0L,
          all(custom_result$summary$primary_test == "combined_fisher"))
for (j in seq_along(pathways)) {
  got <- custom_result$details[[j]]$competitive
  manual <- t(vapply(c(list(match(pathways[[j]], colnames(Y))), got$matched_sets),
                      custom, numeric(2L)))
  close(got$scores, manual)
}

# Numeric pathways refer to original bundle columns. Removing a constant gene
# must not silently shift a numeric pathway to different genes.
with_constant <- bundle
with_constant$Y[, 1L] <- 0
index_statistic <- function(idx) c(index_sum = sum(idx))
numeric_path <- passage_analyze_bundle(with_constant, pathways = list(target = 2:9),
    B = 3L, seed = seed, gene_bins = bins[-1L, ], min_size = 5L, verbose = FALSE,
    statistics = index_statistic)
named_path <- passage_analyze_bundle(with_constant,
    pathways = list(target = colnames(Y)[2:9]), B = 3L, seed = seed,
    gene_bins = bins[-1L, ], min_size = 5L, verbose = FALSE, statistics = index_statistic)
close(numeric_path$details[[1L]]$competitive$scores,
      named_path$details[[1L]]$competitive$scores)
stopifnot(numeric_path$summary$primary_statistic == sum(1:8),
          numeric_path$coverage$n_input == 8L, numeric_path$coverage$n_observed == 8L)

# A matching failure is detected before fitting or calling a custom statistic.
bad_bins <- bins
bad_bins$bin[1L] <- "unsupported"
callback_called <- FALSE
never <- function(idx) { callback_called <<- TRUE; 1 }
expect_error(passage_analyze_bundle(bundle, B = 3L, gene_bins = bad_bins,
    min_size = 5L, verbose = FALSE, statistics = never))
stopifnot(!callback_called)

cat("Bundle integration: legacy default, shared statistic panel, combinations, rank grid, local GP, custom callbacks, covariates and diagnostics passed\n")
