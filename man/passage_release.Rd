\name{passage_analyze_bundle}
\alias{passage_analyze_bundle}
\alias{passage_validate_bundle}
\alias{passage_sample_unique_sets}
\alias{passage_subspace_test}
\title{Validated Bundle Analysis and Spatial Subspace Reference Tests}
\description{Validate portable spatial expression inputs, run competitive pathway analysis with unique matched sets, or run the scale-free trace and largest-root reference tests.}
\usage{
passage_validate_bundle(bundle, min_genes = 5L)
passage_sample_unique_sets(pathway, strata, B = 9999L, seed = 1L,
  reference = c("complement", "universe"))
passage_analyze_bundle(bundle, pathways = NULL, B = 9999L, seed = 1L,
  K = 6L, m = 20L, range_grid = NULL, gene_bins = NULL,
  min_size = 15L, max_size = 500L, verbose = TRUE,
  statistics = "score_z", combine = "none", rank_method = "fixed",
  variance_threshold = 0.90, max_K = NULL, min_K = 2L,
  factor_method = "pca", scale = FALSE, rank_grid = NULL,
  estimation = c("global", "pathway"), local_K = NULL,
  local_variance_threshold = 0.95, local_method = "pca",
  gp_aggregations = c("sum", "max"), gp_top_k = 2L,
  gp_kernel_method = "lowrank", gp_n_features = 64L,
  gp_ranges = c(0.05, 0.1, 0.2, 0.4), keep_matched_sets = FALSE,
  reference = c("complement", "universe"))
passage_subspace_test(Y, coords = NULL, X = NULL, spatial_kernel = NULL,
  kernel_range = 0.25, variance_threshold = 0.95, K = NULL,
  B = 9999L, seed = 1L, max_n = 3000L)
}
\arguments{
\item{bundle}{A list or RDS path with Y, coords, X, aligned row names, and optional cell_type_proportions and pathway_sets.}
\item{min_genes}{Minimum number of genes in a bundle.}
\item{pathway}{Integer indices of pathway genes.}
\item{strata}{One matching stratum per background gene.}
\item{reference}{Legacy target-excluded complement or full-universe conditional gene-label randomization; see passage_competitive_test. Universe allows target overlap while retaining unique genes within each set.}
\item{B}{Positive integer number of Monte Carlo draws.}
\item{seed}{Random seed.}
\item{pathways}{Named list of gene-symbol vectors.}
\item{K}{Factor rank; NULL in the subspace test selects by explained variance.}
\item{m}{Vecchia neighbor count.}
\item{range_grid}{Prespecified candidate spatial ranges.}
\item{gene_bins}{Data frame with gene and bin, optionally module, aligned to nonconstant Y columns.}
\item{min_size, max_size}{Observed pathway-size limits.}
\item{verbose}{Print progress.}
\item{Y}{Locations by genes expression matrix.}
\item{coords}{Aligned spatial coordinates.}
\item{X}{Fixed, full-rank adjustment design. An intercept is included.}
\item{spatial_kernel}{Optional positive semidefinite spatial kernel matrix.}
\item{kernel_range}{Matern 3/2 range in coordinate units.}
\item{variance_threshold}{Fraction of variance retained for global variance-based rank selection or subspace testing.}
\item{max_n}{Safety limit for the dense reference calculation.}
\item{statistics}{Vector of built-in global statistics, or a custom function of gene indices returning a named numeric vector. Default score_z is retained. For pathway estimation use gp_aggregations and leave statistics at its default.}
\item{combine}{none, cauchy, minp, fisher, or a vector of combination choices. Shared-background calibration accounts for component dependence.}
\item{rank_method}{fixed or variance global rank selection.}
\item{max_K, min_K}{Bounds on global rank selection.}
\item{factor_method}{Global factor representation accepted by passage_fit_engine_pca.}
\item{scale}{Whether to standardize genes during global fitting and scoring. This is not a spot-level heteroskedasticity correction.}
\item{rank_grid}{Optional prespecified distinct global fixed ranks. Each statistic/rank pair becomes a component of the shared test.}
\item{estimation}{global fits a shared engine; pathway refits PCA/FA for every target and background set. A custom statistic defines its own estimation.}
\item{local_K, local_variance_threshold, local_method}{Pathway-only fixed rank, PCA variance-retention threshold, and pca/factor_analysis choice.}
\item{gp_aggregations}{sum, max and/or topk combinations of local GP component evidence.}
\item{gp_top_k}{Prespecified component count for topk.}
\item{gp_kernel_method, gp_n_features, gp_ranges}{Local GP cache approximation, feature count and spatial ranges in the internally normalized coordinate units.}
\item{keep_matched_sets}{Retain full sampled gene-index lists for each pathway. FALSE avoids this potentially large storage cost. Seeds, score matrices and matching diagnostics are retained.}
}
\value{Validation returns the checked bundle. Sampling returns a list of integer sets. Analysis returns the primary summary, all component/combination tests in a long tests table, fit diagnostics in estimation, shared scores and matching diagnostics in details, coverage, configuration and session information. The primary test is the first requested combination, otherwise the first statistic; choose this ordering before examining results. Subspace testing returns statistics, Monte Carlo p-values, null draws, rank and assumptions.}
\details{Competitive sets preserve stratum counts and contain no duplicate genes within a set. The default complement reference excludes tested genes; universe allows target overlap under a conditional uniform-gene-label null. Insufficient support is an error. The default matching bins use expression, detection and variance; module matching requires a supplied module column. All tested pathways belong to one BH/BY family. Valid multiplicity correction requires valid marginal p-values; study-specific calibration remains required. Driver ranks and PCA-based variance metrics are exploratory. The subspace reference null requires independent Gaussian rows conditional on a fixed design and permits arbitrary gene covariance. Its two component tests are combined by Bonferroni. It is self-contained, not competitive. Neither API establishes biological validity of supplied cell fractions.}
\seealso{\code{\link{passage_competitive_test}}, \code{\link{passage_make_statistic}}, \code{\link{passage_make_local_gp_statistic}}}
