\name{passage_analyze_bundle}
\alias{passage_analyze_bundle}
\alias{passage_validate_bundle}
\alias{passage_sample_unique_sets}
\alias{passage_subspace_test}
\title{Validated Bundle Analysis and Spatial Subspace Reference Tests}
\description{Validate portable spatial expression inputs, run competitive pathway analysis with unique matched sets, or run the scale-free trace and largest-root reference tests.}
\usage{
passage_validate_bundle(bundle, min_genes = 5L)
passage_sample_unique_sets(pathway, strata, B = 9999L, seed = 1L)
passage_analyze_bundle(bundle, pathways = NULL, B = 9999L, seed = 1L,
  K = 6L, m = 20L, range_grid = NULL, gene_bins = NULL,
  min_size = 15L, max_size = 500L, verbose = TRUE)
passage_subspace_test(Y, coords = NULL, X = NULL, spatial_kernel = NULL,
  kernel_range = 0.25, variance_threshold = 0.95, K = NULL,
  B = 9999L, seed = 1L, max_n = 3000L)
}
\arguments{
\item{bundle}{A list or RDS path with Y, coords, X, aligned row names, and optional cell_type_proportions and pathway_sets.}
\item{min_genes}{Minimum number of genes in a bundle.}
\item{pathway}{Integer indices of pathway genes.}
\item{strata}{One matching stratum per background gene.}
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
\item{variance_threshold}{Fraction of variance retained for subspace testing.}
\item{max_n}{Safety limit for the dense reference calculation.}
}
\value{Validation returns the checked bundle. Sampling returns a list of integer sets. Analysis returns summaries, null draws, coverage, configuration and session information. Subspace testing returns statistics, Monte Carlo p-values, null draws, rank and assumptions.}
\details{Competitive sets preserve stratum counts and contain no duplicates or tested genes. Insufficient support is an error. The default matching bins use expression, detection and variance; module matching requires a supplied module column. All tested pathways belong to one BH/BY family. Valid multiplicity correction requires valid marginal p-values; study-specific calibration remains required. Driver ranks and PCA-based variance metrics are exploratory. The subspace reference null requires independent Gaussian rows conditional on a fixed design and permits arbitrary gene covariance. Its two component tests are combined by Bonferroni. It is self-contained, not competitive. Neither API establishes biological validity of supplied cell fractions.}
