if (dir.exists("R")) {
  for (f in sort(list.files(file.path("R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

set.seed(20260524)

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 10),
  y = seq(0, 1, length.out = 10)
))
n <- nrow(coords)
g <- 60
Y <- matrix(rnorm(n * g, sd = 0.8), nrow = n, ncol = g)
colnames(Y) <- paste0("g", seq_len(g))

d <- as.matrix(dist(coords))
K <- exp(-d / 0.25) + diag(1e-6, n)
v <- as.numeric(t(chol(K)) %*% rnorm(n))
u <- rnorm(n)
for (jj in 1:8) {
  Y[, jj] <- Y[, jj] + 1.2 * v + 0.6 * u
}
for (jj in 9:14) {
  Y[, jj] <- Y[, jj] + 0.6 * u
}

engine <- passage_fit_engine_pca(
  Y = Y,
  coords = coords,
  K = 3,
  m = 8,
  range_grid = c(0.10, 0.20, 0.35),
  verbose = FALSE
)
pre <- passage_h_precompute(engine, X = matrix(1, n, 1))
bins <- passage_make_gene_bins(Y, n_mean_bins = 4, n_detect_bins = 4)
bins$module <- paste0("m", rep(1:6, length.out = nrow(bins)))

signal <- paste0("g", 1:10)
null <- paste0("g", 31:45)

m_signal <- passage_pathway_covariance_metrics(engine, Y, signal, precomp = pre)
m_null <- passage_pathway_covariance_metrics(engine, Y, null, precomp = pre)
stopifnot(is.finite(m_signal$summary[["cwPVE_trace"]]))
stopifnot(is.finite(m_signal$summary[["cEPSV"]]))
stopifnot(is.finite(m_signal$summary[["pc1_spatial_fraction"]]))
stopifnot(m_signal$summary[["mean_propSV_conditional"]] > m_null$summary[["mean_propSV_conditional"]])

gene_stats <- passage_gene_spatial_stats(engine, Y, precomp = pre)
z <- passage_competitive_gene_stat_z(signal, gene_stats$propSV, bins, colnames(Y))
stopifnot(is.finite(z$statistic))
stopifnot(is.finite(z$p))

cc <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "mean_propSV_conditional", B = 25, seed = 11
)
stopifnot(is.finite(cc$analytic$p))
stopifnot(is.finite(cc$permutation$p))

ce <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "cEPSV", B = 10, seed = 12
)
stopifnot(is.finite(ce$permutation$p))

sq <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "score_Q", B = 10, seed = 112
)
stopifnot(is.finite(sq$permutation$p))
sz <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "score_z", B = 10, seed = 113
)
stopifnot(is.finite(sz$permutation$p))
stopifnot(is.finite(passage_pathway_score_stat(engine, Y, signal, pre, output = "z")))
gz <- passage_gene_score_z(engine, Y, pre, gene_names = colnames(Y))
stopifnot(nrow(gz) == ncol(Y))
stopifnot(is.finite(gz$score_z_gene[[1]]))
rz <- passage_robust_pathway_score_stat(engine, Y, signal, pre, gene_names = colnames(Y), gene_scores = gz)
stopifnot(is.finite(rz))
residualize_with_qr <- if (exists("passage_residualize_with_qr", mode = "function")) {
  passage_residualize_with_qr
} else {
  getFromNamespace("passage_residualize_with_qr", "PASSAGE")
}
R_all <- residualize_with_qr(Y, pre$qr)
rz_fast <- passage_robust_pathway_score_stat(
  engine, Y, signal, pre, gene_names = colnames(Y), gene_scores = gz,
  residual_matrix = R_all
)
stopifnot(abs(rz_fast - rz) < 1e-10)
sr <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "score_robust_z", B = 10, seed = 114
)
stopifnot(is.finite(sr$permutation$p))
pa <- passage_gene_pERSA(engine, Y, pre, gene_names = colnames(Y), gene_bins = bins)
stopifnot(nrow(pa) == ncol(Y))
stopifnot(is.finite(pa$pERSA_z[[1]]))
pa_rank <- passage_gene_pERSA(
  engine, Y, pre, gene_names = colnames(Y), gene_bins = bins,
  standardization = "rank"
)
stopifnot(is.finite(pa_rank$pERSA_z[[1]]))
pa_active <- passage_gene_pERSA(
  engine, Y, pre, gene_names = colnames(Y), gene_bins = bins,
  standardization = "active"
)
stopifnot(is.finite(pa_active$pERSA_z[[1]]))
pz <- passage_pERSA_pathway_score_stat(
  engine, Y, signal, pre, gene_names = colnames(Y), gene_activity = pa,
  gene_bins = bins, residual_matrix = R_all
)
stopifnot(is.finite(pz))
sp <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "score_pERSA_z", B = 10, seed = 115
)
stopifnot(is.finite(sp$permutation$p))
for (st in c("score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z")) {
  tmp <- passage_conditional_competitive_test(
    engine, Y, signal, gene_bins = bins, precomp = pre,
    statistic = st, B = 6, seed = 120 + match(st, c("score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z"))
  )
  stopifnot(is.finite(tmp$permutation$p))
}
sp_topk <- passage_sparse_topk_pathway_score_stat(
  engine, Y, signal, pre, gene_names = colnames(Y), gene_scores = gz,
  residual_matrix = R_all
)
stopifnot(is.finite(sp_topk))
lse <- passage_gene_lse(gz$score_z_gene, gene_bins = bins)
stopifnot(is.finite(lse$lse_z[[1]]))
lse_stat <- passage_gene_activity_pathway_score_stat(
  lse$lse_z, signal, Y, pre, gene_names = colnames(Y),
  residual_matrix = R_all
)
stopifnot(is.finite(lse_stat))
hsic_stat <- passage_pathway_factor_hsic_stat(
  engine, Y, signal, pre, gene_names = colnames(Y),
  residual_matrix = R_all
)
stopifnot(is.finite(hsic_stat))
hot_stat <- passage_pathway_activity_hotspot_stat(
  engine, Y, signal, pre, gene_names = colnames(Y),
  residual_matrix = R_all
)
stopifnot(is.finite(hot_stat))
for (st in c("score_sparse_topk_z", "score_lse_z", "score_factor_hsic",
             "score_coherence", "score_activity_hotspot")) {
  tmp <- passage_conditional_competitive_test(
    engine, Y, signal, gene_bins = bins, precomp = pre,
    statistic = st, B = 6, seed = 130 + match(st, c("score_sparse_topk_z", "score_lse_z", "score_factor_hsic",
                                                    "score_coherence", "score_activity_hotspot"))
  )
  stopifnot(is.finite(tmp$permutation$p))
}

for (sampler in c("independent", "module", "mcmc", "importance", "knockoff")) {
  ss <- passage_conditional_competitive_test(
    engine, Y, signal, gene_bins = bins, precomp = pre,
    statistic = "score_z", B = 8, seed = 210 + match(sampler, c("independent", "module", "mcmc", "importance", "knockoff")),
    sampler = sampler, mcmc_burnin = 5, mcmc_thin = 2,
    importance_oversample = 3
  )
  stopifnot(is.finite(ss$permutation$p))
  stopifnot(identical(ss$permutation$sampler, sampler))
}

balance_fun <- passage_make_expression_coherence_balance(
  Y, precomp = pre, metrics = c("expr_mean_abs_cor", "expr_pc1_fraction")
)
balance_target <- balance_fun(passage_resolve_pathway(signal, colnames(Y)))
balance_tol <- pmax(c(expr_mean_abs_cor = 0.25, expr_pc1_fraction = 0.25) * balance_target, 0.05)
ce_bal <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "cEPSV", B = 10, seed = 121,
  balance_fun = balance_fun,
  balance_target = balance_target,
  balance_tol = balance_tol,
  balance_max_tries = 10
)
stopifnot(is.finite(ce_bal$permutation$p))
stopifnot(!is.null(ce_bal$permutation$balance))
stopifnot(is.finite(ce_bal$permutation$balance$acceptance_rate))

fast_ctx <- passage_competitive_fast_context(engine, Y, precomp = pre)
factor_balance_fun <- passage_make_factor_coherence_balance(fast_ctx)
factor_balance_target <- factor_balance_fun(passage_resolve_pathway(signal, colnames(Y)))
factor_balance_tol <- pmax(c(loading_pc1_fraction = 0.25,
                             loading_mean_factor_coherence = 0.25) * factor_balance_target, 0.05)
ce_factor_bal <- passage_conditional_competitive_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "cEPSV", B = 10, seed = 122,
  balance_fun = factor_balance_fun,
  balance_target = factor_balance_target,
  balance_tol = factor_balance_tol,
  balance_max_tries = 10
)
stopifnot(is.finite(ce_factor_bal$permutation$p))
stopifnot(!is.null(ce_factor_bal$permutation$balance))

coh <- passage_pathway_coherence_test(
  engine, Y, signal, gene_bins = bins, precomp = pre,
  statistic = "pc1_spatial_fraction", B = 10, seed = 13
)
stopifnot(is.finite(coh$permutation$p))

region <- factor(ifelse(coords[, 1] > stats::median(coords[, 1]), "right", "left"),
                 levels = c("left", "right"))
re <- passage_region_enrichment_test(
  engine, Y, signal, region = region, precomp = pre,
  n_perm = 10, seed = 14
)
stopifnot(is.finite(re$p))
stopifnot(ncol(re$Q_region) == 2L)

cat("PASSAGE competitive metrics test passed\n")
