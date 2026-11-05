if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}

expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
close <- function(a, b, tol = 1e-10) {
  stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE)))
}

# A user statistic is the only part that knows what expression score means.
# The testing layer must reuse it without changing the chosen pathway rule.
genes <- paste0("custom_gene_", seq_len(12L))
value <- c(4, 7, 1, 2, 4, 8, 3, 9, 2, 6, 5, 10)
sets <- list(c(3L, 4L), c(5L, 6L), c(7L, 8L), c(9L, 10L), c(11L, 12L))
calls <- list()
statistic <- function(idx) {
  calls[[length(calls) + 1L]] <<- idx
  c(total = sum(value[idx]), peak = max(value[idx]))
}
ans <- passage_competitive_test(genes[1:2], statistic, gene_names = genes,
  matched_sets = sets, combine = c("cauchy", "minp", "fisher"))
expected <- t(vapply(c(list(1:2), sets), function(idx)
  c(total = sum(value[idx]), peak = max(value[idx])), numeric(2)))
close(ans$scores, expected)
stopifnot(length(calls) == 6L, identical(as.integer(calls[[1L]]), 1:2),
  identical(ans$components$statistic, colnames(expected)),
  all(ans$components$B == 5L), all(is.na(ans$components$mc_lower)),
  all(is.na(ans$components$mc_upper)),
  setequal(ans$combinations$method, c("cauchy", "minp", "fisher")))
exceed <- colSums(sweep(expected[-1L, , drop = FALSE], 2L, expected[1L, ], ">="))
close(ans$components$p_competitive, (1 + exceed) / nrow(expected))
close(ans$components$exceedances, exceed)
close(ans$components$observed, expected[1L, ])

# Independently form pooled upper-tail probabilities and combine them. The
# final combined result is an empirical tail, never an independence chi-square
# reference or an uncalibrated minimum p-value.
pooled <- apply(expected, 2L, function(x) vapply(x, function(v) mean(x >= v), numeric(1)))
manual <- list(minp = -apply(pooled, 1L, min),
               fisher = -2 * rowSums(log(pooled)))
clipped <- pmin(pmax(pooled, 0.5 / nrow(expected)), 1 - 0.5 / nrow(expected))
manual$cauchy <- rowMeans(tan(pi * (0.5 - clipped)))
for (method in names(manual)) {
  z <- manual[[method]]
  got <- ans$combinations$p_competitive[ans$combinations$method == method]
  close(got, mean(z >= z[1L]))
}
stopifnot(all(is.na(ans$combinations$mc_lower)), all(is.na(ans$combinations$mc_upper)))

# Lower-tail orientation must apply before component ranking and combination.
lower <- passage_competitive_test(1:2, function(idx) c(total = -sum(value[idx]), peak = max(value[idx])),
  gene_names = genes, matched_sets = sets,
  alternative = c(peak = "greater", total = "less"), combine = c("cauchy", "minp", "fisher"))
close(lower$components$p_competitive, ans$components$p_competitive)
close(lower$combinations$p_competitive, ans$combinations$p_competitive)

# Ties count toward the tail, including statistics that are entirely constant.
ties <- passage_competitive_test(1:2, function(idx) c(constant = 3, another = -5),
  gene_names = genes, matched_sets = sets, combine = c("cauchy", "minp", "fisher"))
stopifnot(all(ties$components$p_competitive == 1), all(ties$combinations$p_competitive == 1))
duplicated_stat <- passage_competitive_test(1:2,
  function(idx) c(first = sum(value[idx]), duplicate = sum(value[idx])),
  gene_names = genes, matched_sets = sets, combine = c("cauchy", "minp", "fisher"))
close(duplicated_stat$combinations$p_competitive,
      rep(duplicated_stat$components$p_competitive[1L], 3L))

# Exhaustive designation of each row as the target checks the exchangeable-row
# rank law without noisy simulation: at most floor(alpha*N) targets reject.
# This is a test of calibration algebra, not proof of gene-matching validity.
fixed_scores <- cbind(first = (1:20)^2, second = c(11:20, 1:10),
                      tied = rep(1:5, each = 4L))
all_p <- vapply(seq_len(nrow(fixed_scores)), function(i) {
  z <- passage_calibrate_competitive_scores(fixed_scores[c(i, setdiff(1:20, i)), , drop = FALSE],
    combine = c("cauchy", "minp", "fisher"), mc_intervals = FALSE)
  c(z$components$p_competitive, z$combinations$p_competitive)
}, numeric(6))
for (alpha in c(.05, .10, .25, .50))
  stopifnot(all(rowSums(all_p <= alpha) <= floor(alpha * 20L + 1e-8)))
permutation <- c(1L, 6L, 3L, 2L, 5L, 4L)
reordered <- passage_calibrate_competitive_scores(expected[permutation, , drop = FALSE],
  combine = c("cauchy", "minp", "fisher"), mc_intervals = FALSE)
close(reordered$components$p_competitive, ans$components$p_competitive)
close(reordered$combinations$p_competitive, ans$combinations$p_competitive)

# Generated sets preserve bin AND module margins, exclude the target, contain
# no duplicated gene, and leave the caller's RNG stream untouched.
gene_names <- paste0("stratum_gene_", seq_len(80L))
bins <- data.frame(gene = gene_names, bin = rep(c("A", "B"), each = 40L),
                   module = rep(rep(c("M1", "M2"), each = 20L), 2L))
target <- c(1L, 21L, 41L, 61L)
seen <- list()
record_sets <- function(idx) {
  seen[[length(seen) + 1L]] <<- idx
  c(total = sum(idx))
}
set.seed(9346L)
saved_seed <- .Random.seed
sampled <- passage_competitive_test(target, record_sets, gene_bins = bins,
  gene_names = gene_names, B = 39L, seed = 94L)
stopifnot(identical(saved_seed, .Random.seed), length(seen) == 40L,
  is.finite(sampled$components$mc_lower), is.finite(sampled$components$mc_upper),
  sampled$components$mc_lower <= sampled$components$mc_upper)
strata <- paste(bins$bin, bins$module, sep = "|")
target_margin <- table(factor(strata[target], levels = unique(strata)))
stopifnot(all(vapply(seen[-1L], function(idx) length(idx) == 4L && !anyDuplicated(idx) &&
  !any(idx %in% target) && identical(table(factor(strata[idx], levels = unique(strata))), target_margin), logical(1))))
repeat_sample <- passage_competitive_test(target, function(idx) c(total = sum(idx)),
  gene_bins = bins, gene_names = gene_names, B = 39L, seed = 94L)
close(sampled$scores, repeat_sample$scores)
rm(".Random.seed", envir = .GlobalEnv)
invisible(passage_competitive_test(target, function(idx) c(total = sum(idx)),
  gene_bins = bins, gene_names = gene_names, B = 5L, seed = 95L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assign(".Random.seed", saved_seed, envir = .GlobalEnv)

# Insufficient support and malformed scores must fail closed, without quietly
# discarding a component, changing strata, or inventing a finite replacement.
bad_bins <- bins
bad_bins$bin[1L] <- "unique"
expect_error(passage_competitive_test(target, function(idx) c(total = sum(idx)),
  gene_names = gene_names, gene_bins = bad_bins, B = 5L))
expect_error(passage_competitive_test(1:2, function(idx) c(bad = NA_real_),
  gene_names = genes, matched_sets = sets))
expect_error(passage_competitive_test(1:2, function(idx) c(bad = Inf),
  gene_names = genes, matched_sets = sets))
expect_error(passage_competitive_test(1:2, function(idx) if (1L %in% idx) c(a = 1) else c(b = 1),
  gene_names = genes, matched_sets = sets))
expect_error(passage_competitive_test(1:2, function(idx) c(a = 1), gene_names = genes,
  matched_sets = list(c(1L, 3L))))
expect_error(passage_competitive_test(1:2, function(idx) c(a = 1), gene_names = genes,
  matched_sets = list(c(3L, 3L))))
expect_error(passage_competitive_test(1:2, function(idx) c(a = 1), gene_names = genes,
  B = 0L))
expect_error(passage_calibrate_competitive_scores(expected, combine = "unknown"))
expect_error(passage_calibrate_competitive_scores(expected, alternative = "unknown"))

# The factory is only a score calculator. Verify it against the established
# definitions, then run each available statistic through the generic tester.
set.seed(9350L)
coords <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6L), y = seq(0, 1, length.out = 6L)))
X <- cbind(intercept = 1, nuisance = rnorm(nrow(coords)))
Y <- matrix(rnorm(nrow(coords) * 90L), nrow(coords), 90L) +
  X[, 2L] %o% rnorm(90L)
colnames(Y) <- paste0("fit_gene_", seq_len(ncol(Y)))
engine <- passage_fit_engine_pca(Y, coords, X = X, K = 3L, m = 8L,
  range_grid = c(.1, .3), verbose = FALSE)
pre <- passage_h_precompute(engine, X = X)
four <- c("score_z", "score_activity_hotspot", "score_robust_z", "score_sparse_topk_z")
factory <- passage_make_statistic(engine, Y, statistics = four, precomp = pre)
observed <- factory(1:10)
reference <- c(
  score_z = passage_pathway_score_stat(engine, Y, 1:10, pre, output = "z"),
  score_activity_hotspot = passage_pathway_activity_hotspot_stat(engine, Y, 1:10, pre),
  score_robust_z = passage_robust_pathway_score_stat(engine, Y, 1:10, pre),
  score_sparse_topk_z = passage_sparse_topk_pathway_score_stat(engine, Y, 1:10, pre))
close(observed, reference, 1e-8)
all_statistics <- c(four, "score_lse_z", "score_factor_hsic", "score_coherence",
  "score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z",
  "score_Q", "mean_propSV_conditional", "cEPSV", "ePSV", "cwPVE_trace", "cwPVE_top")
all_factory <- passage_make_statistic(engine, Y, statistics = all_statistics, precomp = pre)
all_result <- passage_competitive_test(1:10, all_factory, gene_names = colnames(Y),
  B = 3L, seed = 95L, combine = c("cauchy", "minp", "fisher"))
stopifnot(ncol(all_result$scores) == 17L, all(is.finite(all_result$scores)),
  identical(all_result$components$statistic, all_statistics))

cat("Competitive interface: generic callback, pooled combination calibration, ties, tails, matching, RNG, failure handling and all statistic adapters passed\n")
