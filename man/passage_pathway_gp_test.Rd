\name{passage_test_pathway_gp}
\alias{passage_test_pathway_gp}
\alias{passage_prepare_gp_kernel}
\title{Adaptive Spatial GP Model Comparison for Pathway Components}
\description{Fit spatial covariance parameters to unsmoothed pathway components and test the improvement over nonspatial models. A joint Gaussian rotation reference repeats the parameter search for every null draw. Coordinate-only kernel caches can be reused across pathways with the same locations and adjustment design.}
\usage{
passage_prepare_gp_kernel(coords, X = NULL,
  ranges = c(0.05, 0.1, 0.2, 0.4),
  kernel_method = c("lowrank", "dense"), n_features = 64L,
  max_dense_n = 1000L)
passage_test_pathway_gp(fit, coords = NULL, kernel_cache = NULL,
  ranges = c(0.05, 0.1, 0.2, 0.4),
  fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
  B = 999L, seed = 1L, kernel_method = c("lowrank", "dense"),
  n_features = 64L, max_dense_n = 1000L, batch_size = 64L,
  keep_null = TRUE)
}
\arguments{
\item{coords}{Finite aligned spatial coordinates. Coordinate units are retained. For unit-domain ranges, translate the coordinates and divide all axes by the same maximum axis span. Supplied spot IDs must match the fitted scores exactly.}
\item{X}{Fixed full-rank adjustment design; an intercept is included. A reused cache must have exactly the same adjustment design and spot order as the fit.}
\item{ranges}{Prespecified positive Mat\'ern 3/2 range grid in coordinate units. The complete grid is searched for the observed scores and every null draw.}
\item{kernel_method}{Coordinate-only low-rank Nystr\"om approximation, or a dense reference. Both construct positive-semidefinite covariance kernels.}
\item{n_features}{Maximum number of coordinate-selected landmarks for low-rank kernels. Increasing it can improve covariance approximation at higher computational cost.}
\item{max_dense_n}{Spot limit for the dense reference. No automatic subsampling occurs.}
\item{fit}{An unsmoothed \code{passage_fit_pathway()} result. Smoothed scores and descriptive GpGp result objects are not inputs.}
\item{kernel_cache}{Optional \code{passage_prepare_gp_kernel()} result. Cache settings govern the calculation; explicitly conflicting settings are rejected.}
\item{fractions}{Prespecified spatial variance fraction grid in \eqn{[0,1)}, including zero for the null. Fractions refer to average residual-space covariance, after kernel trace normalization.}
\item{B}{Number of null rotations. The minimum Monte Carlo p-value is \eqn{1/(B+1)}.}
\item{seed}{Random seed; the prior RNG state is restored.}
\item{batch_size}{Number of null draws processed per batch. Smaller batches reduce working memory.}
\item{keep_null}{Whether to retain the scalar pathway null statistics in the result.}
}
\value{
\code{passage_prepare_gp_kernel()} returns a coordinate/design-specific cache with projected kernel spectra and settings.

\code{passage_test_pathway_gp()} returns the scalar sum of twice the profile restricted-loglikelihood improvements in \code{statistic}, a Monte Carlo \code{p}, and \code{component_fits}. Each component row reports its statistic, fitted grid range and fraction, total/spatial/noise variance estimates, and grid-boundary flags. The result also contains method/calibration labels, optional null statistics, Monte Carlo count and seed, residual degrees of freedom, retained rank, configuration and timing.
}
\details{
After adjustment, each component is modeled with covariance proportional to
\eqn{(1-h)I + hK_r}, where \eqn{h} is the spatial variance fraction and
\eqn{K_r} is the projected Mat\'ern covariance at range \eqn{r}, normalized to have trace equal to the residual degrees of freedom. The calculation profiles out total component variance, maximizes over the supplied fraction and range grids, and sums the marginal component improvements over \eqn{h=0}. For correlated FA scores this sum is a composite statistic rather than an independent-factor joint likelihood ratio. No chi-square reference is used.

The null jointly rotates the unsmoothed score matrix within the covariate residual space, preserving the complete component Gram matrix. Every draw repeats the same range/fraction maximization. Under Gaussian row errors with identity row covariance conditional on a fixed design, pathway PCA and the implemented FA estimator are equivariant to these row rotations; expression Gram matrices, loadings and selected ranks can be conditioned on without numerically refitting them. Coordinate-selected landmarks and the prespecified parameter grids remain fixed across the observed and null calculations.

The low-rank option avoids dense spot covariance decompositions by using Nystr\"om features and projected small-matrix spectra. It approximates the spatial covariance model and can lose power for patterns outside the retained features. Use feature-count sensitivity analyses and the dense reference where feasible. Kernel construction may be reused only across fits with aligned spots and the same adjustment design; never choose grids or landmarks using an unaccounted expression-dependent search.

Grid parameter estimates are descriptive plug-in estimates without confidence intervals. A zero fraction gives no identified spatial range; upper-fraction and range-boundary flags indicate where a wider prespecified grid may be useful in a future sensitivity analysis. They do not justify selecting a new primary analysis after seeing its p-value.

The test is self-contained for the retained pathway representation. It does not implement matched-set competitive testing, PC-wise SKAT, ordinary Tracy--Widom testing, or across-pathway multiplicity correction. Normalized counts, heteroskedastic spots and expression-estimated covariates require study-specific calibration checks. GpGp is not required by this testing option; \code{passage_fit_pathway_gp()} remains a separate descriptive fitting API.
}
\seealso{\code{\link{passage_fit_pathway}}, \code{\link{passage_test_pathway_fit}}, \code{\link{passage_fit_pathway_gp}}, \code{\link{passage_analyze_bundle}}}
\examples{
set.seed(7)
Y <- matrix(rnorm(50 * 12), 50, 12,
  dimnames = list(paste0("s", 1:50), paste0("g", 1:12)))
xy <- cbind(runif(50), runif(50)); rownames(xy) <- rownames(Y)
fit <- passage_fit_pathway(Y, paste0("g", 1:8), K = 2)
cache <- passage_prepare_gp_kernel(xy, X = fit$X,
  ranges = c(0.1, 0.3), n_features = 12)
ans <- passage_test_pathway_gp(fit, kernel_cache = cache, B = 19)
ans$p
ans$component_fits
}
