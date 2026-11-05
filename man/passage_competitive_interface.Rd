\name{passage_competitive_test}
\alias{passage_competitive_test}
\alias{passage_matching_diagnostics}
\alias{passage_calibrate_competitive_scores}
\alias{passage_make_statistic}
\alias{passage_make_local_gp_statistic}
\title{Statistic-Independent Matched-Set Competitive Tests}
\description{Construct reusable global or pathway-local statistic functions, inspect matching support, and calibrate arbitrary numeric pathway statistics against shared matched gene sets.}
\usage{
passage_matching_diagnostics(pathway, gene_bins, gene_names = NULL,
  reference = c("complement", "universe"))
passage_calibrate_competitive_scores(scores, combine = "none",
  alternative = "greater", mc_intervals = FALSE, conf_level = 0.95)
passage_competitive_test(pathway, statistic, gene_bins = NULL,
  gene_names = NULL, B = 9999L, seed = 1L, matched_sets = NULL,
  combine = "none", alternative = "greater", conf_level = 0.95,
  keep_matched_sets = TRUE, reference = c("complement", "universe"))
passage_make_statistic(engine, Y, statistics = "score_z",
  precomp = NULL, X = NULL, gene_bins = NULL,
  score_weight_scheme = "equal", covariance = "shrink",
  coherence = "signed", gene_names = NULL, ...)
passage_make_local_gp_statistic(Y, coords, X = NULL, method = "pca",
  K = NULL, variance_threshold = 0.95, kernel_cache = NULL,
  ranges = c(0.05, 0.1, 0.2, 0.4), kernel_method = "lowrank",
  n_features = 64L, aggregations = c("sum", "max"), top_k = 2L,
  fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
  seed = 1L)
}
\arguments{
\item{pathway}{Gene names or integer indices in the declared background universe.}
\item{gene_bins}{Data frame with aligned gene and bin columns and optionally module. Without bins, all supplied genes form one matching stratum.}
\item{gene_names}{Unique background gene names, in expression-column order.}
\item{scores}{Finite numeric matrix, with observed pathway in row one, shared backgrounds in later rows, and named statistics in columns.}
\item{combine}{One or more of cauchy, minp, fisher; or none. Each combination is separately calibrated against the shared background rows.}
\item{alternative}{greater or less, either one value or one per statistic; named vectors must match the statistic names.}
\item{mc_intervals}{Whether rows represent independent background draws conditional on a fixed target, permitting binomial intervals for individual statistics. The generic test sets this automatically for internally drawn sets only.}
\item{conf_level}{Confidence level for individual conditional background-tail probabilities.}
\item{statistic}{Function accepting a gene-index vector and returning one number or a uniquely named numeric vector. It must use identical rules for target and background sets, including any pathway-specific fitting or selection.}
\item{B}{Positive integer number of background draws. Supplied matched_sets determines the count when present.}
\item{seed}{Random seed. Generic competitive inference preserves the caller's random state; a NULL seed uses the current stream. The local fitting factory requires a nonnegative integer seed.}
\item{matched_sets}{Optional gene-index lists with unique genes and the declared stratum composition. Complement references must exclude target genes; universe references may overlap the target or equal it.}
\item{reference}{Complement retains the legacy target-excluded comparison. Universe samples uniformly within strata from all eligible genes, including target genes: a conditional gene-label randomization reference. The choice changes the null and must be prespecified.}
\item{keep_matched_sets}{Retain sampled gene-index lists. FALSE reduces memory for large pathway collections; score matrices and matching summaries remain available.}
\item{engine}{Fitted global PASSAGE engine.}
\item{Y}{Finite spots-by-genes expression matrix aligned to the engine.}
\item{statistics}{Unique built-in statistic names. All 17 choices in passage_conditional_competitive_test are supported, plus mean_abs_spatial_cor.}
\item{precomp}{Optional kernel moments and design from passage_h_precompute.}
\item{X}{Fixed adjustment design. Global defaults inherit the fitting engine's design.}
\item{score_weight_scheme}{equal, var, range or invvar factor weights.}
\item{covariance}{shrink, sample or diagonal residual covariance estimate.}
\item{coherence}{signed, absolute or none loading-coherence weighting.}
\item{...}{Named statistic options: topk_grid, winsor_q, corr, n_pc, scale_floor_frac, shrinkage, ridge, active_tail, ridge_frac. These are fixed before testing.}
\item{coords}{Spatial coordinates aligned to expression rows.}
\item{method}{pca or factor_analysis, fitted separately within every set.}
\item{K}{Fixed local rank, or NULL for variance-based PCA rank. Factor analysis requires K.}
\item{variance_threshold}{PCA variance-retention threshold.}
\item{kernel_cache}{Prepared GP kernel cache with exactly matching coordinates/design.}
\item{ranges}{Prespecified spatial range candidates.}
\item{kernel_method}{lowrank or dense GP approximation.}
\item{n_features}{Number of low-rank features.}
\item{aggregations}{Unique sum, max and/or topk component combinations.}
\item{top_k}{Prespecified number of strongest components summed for topk, capped at the fitted rank.}
\item{fractions}{Prespecified GP spatial-variance fraction candidates, including zero.}
}
\value{The statistic factories return functions of gene indices. Matching diagnostics returns stratum support and matching configuration. Calibration returns component and combination tables, raw scores, pooled ranks, combined scores and calibration metadata. The generic test additionally returns matching diagnostics, sampled sets, configuration and optional local fit diagnostics.}
\details{
For every individual statistic the upper-tail p-value is (1 + exceedances)/(B + 1), after reversing a lower-tail statistic if requested. Exact ties count against the target. Global fitted loadings and GP kernels can be cached; local PCA/FA and component selection are repeated for every gene set. Built-in global scaling follows the units used during estimation.

Combinations use upper-tail ranks of every row among all B+1 rows, including the observed row, then calculate the same combining function for each row and obtain its empirical tail. Cauchy ranks are clipped to [0.5/(B+1), 1-0.5/(B+1)] before transformation to make endpoints finite. Fisher uses minus twice the sum of log ranks; minimum-p uses minus the smallest rank. Neither independence-based chi-square Fisher p-values nor analytic Cauchy-tail approximations are used. Dependence between component statistics is retained. Several returned combinations are separate tests, not permission to select their smallest p-value.

The universe option is a conditional gene-label randomization test: conditional on expression, strata, design and any shared fit, the target must be uniformly distributed among all sets with its stratum counts. Under that assumption, target and independently generated reference sets have exchangeable scores, even with heterogeneous spot noise or strong nuisance structure. Overlap with the target is intentional. No genes repeat within a set; sets can recur across draws. A stratum containing only target genes contributes no randomization. Diagnostics report log reference-set count and whether any variation is possible.

The uniform-label assumption is not guaranteed for curated correlated biological modules. It is distinct from a fixed-pathway, repeated-expression sampling null. Complement exclusion does not supply the same finite-universe randomization argument. Both choices require interpretation under their stated comparison. Module matching changes the comparison question. This option does not correct the self-contained Gaussian rotation reference or guarantee false-discovery control for arbitrary pathway collections.

Binomial Monte Carlo intervals concern individual conditional background-tail probabilities when sets are generated independently. They are not biological confidence intervals and are not supplied for pooled-rank combinations or externally supplied sets. No automatic extra sampling or relaxed matching is performed after seeing a p-value.
}
\examples{
genes <- paste0("g", 1:100)
bins <- data.frame(gene = genes, bin = rep(1:2, each = 50))
gene_measure <- sin(seq_len(100))
score <- function(idx) c(average = mean(gene_measure[idx]),
                         strongest = max(gene_measure[idx]))
ans <- passage_competitive_test(c(1:5, 51:55), score, bins,
  B = 39, seed = 1, combine = "minp")
ans$components
ans$combinations
}
