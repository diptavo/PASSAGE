\name{passage_fit_pathway}
\alias{passage_fit_pathway}
\alias{passage_test_pathway_fit}
\alias{passage_fit_pathway_gp}
\title{Pathway-Only Estimation, Spatial Subspace Testing and Optional GP Fits}
\description{Estimate unsmoothed components using only the selected pathway, test their spatial subspace, and optionally fit descriptive component-wise Gaussian processes. No global PASSAGE engine is used.}
\usage{
passage_fit_pathway(Y, pathway, X = NULL,
  method = c("pca", "factor_analysis"), K = NULL,
  variance_threshold = 0.95, seed = 1L)
passage_test_pathway_fit(fit, coords = NULL, spatial_kernel = NULL,
  kernel_range = 0.25, B = 9999L, seed = 1L, max_n = 3000L)
passage_fit_pathway_gp(fit, coords,
  covfun_name = "matern15_isotropic", m_seq = c(10L, 30L),
  max_iter = 40L, seed = 1L)
}
\arguments{
\item{Y}{Finite spots-by-genes expression matrix with unique gene names. Unnamed spots receive generated IDs.}
\item{pathway}{Gene names or integer column indices. Missing named genes and genes constant after adjustment are reported.}
\item{X}{Fixed, full-rank adjustment design. An intercept is included. If row names are supplied they must exactly match Y.}
\item{method}{Covariance PCA without gene scaling, or maximum-likelihood factor analysis on the gene correlation matrix using unrotated regression scores.}
\item{K}{Prespecified retained rank. NULL selects PCA rank by explained variance. FA requires an explicit identifiable rank and a nonsingular gene correlation matrix. Rank must be below residual degrees of freedom.}
\item{variance_threshold}{Cumulative pathway PCA variance threshold, used only when K is NULL.}
\item{seed}{Random seed; the previous RNG state is restored.}
\item{fit}{An unsmoothed passage_fit_pathway result. A descriptive GP result is not a testing input.}
\item{coords}{Aligned spatial coordinates. Coordinate units are retained. For unit-domain ranges, translate coordinates and divide every axis by the same maximum axis span. GP fitting requires row names identical to fit scores.}
\item{spatial_kernel}{Optional fixed positive-semidefinite kernel. If named, both dimensions must exactly match score IDs.}
\item{kernel_range}{Matern 3/2 range in supplied coordinate units; default 0.25 assumes unit-domain scaling.}
\item{B}{Number of Gaussian-subspace null draws.}
\item{max_n}{Dense-reference spot limit. No silent subsampling is performed.}
\item{covfun_name}{Optional GpGp covariance function.}
\item{m_seq}{Nondecreasing GpGp neighbor counts, capped at the number of spots minus one.}
\item{max_iter}{Maximum GpGp fitting iterations.}
}
\value{Estimation returns a passage_pathway_fit with scores, loadings, rank, covariates, gene coverage, residual variances and diagnostics. Testing returns scale-free trace/root statistics, component and Bonferroni-omnibus Monte Carlo p-values, null draws and assumptions. GP fitting returns separate descriptive fits, parameter estimates, convergence flags and warnings, without p-values or predictions.}
\details{
The fixed-kernel test consumes the estimated unsmoothed score subspace and retains all supplied components. PCA selection is not repeated on a different gene universe. Under iid Gaussian row errors conditional on a fixed design, covariance-only PCA/FA estimation gives a rotation-equivariant score subspace; the Gaussian residual-space reference applies conditional on the selected rank. This does not establish validity for count data, heteroskedasticity, expression-estimated covariates, or data-adaptive spatial kernels. FA regression scores need not be orthogonal: the testing calculation orthonormalizes their span.

GpGp fits are optional and use the same unsmoothed scores. Their covariance estimates are descriptive plug-in estimates; parameter and PCA/FA uncertainty is not propagated. Smoothing before testing introduces spatial row covariance even under a nonspatial null. Using same-data fitted GP values or ranges therefore requires a different calibration that repeats the relevant fitting steps. The fixed-kernel test never uses the optional GP output. For an integrated fitted-covariance test, use \code{\link{passage_test_pathway_gp}}: it fits its own spatial covariance grids and recalibrates that parameter search. It does not reuse GpGp estimates or smoothed predictions.

For the default GpGp Matern 3/2 covariance, GpGp alpha equals PASSAGE range divided by sqrt(3). Actual nugget variance is the fitted process variance times the relative nugget parameter. Inspect all convergence flags, boundary estimates and warnings.

The APIs test one pathway and do not perform across-pathway multiplicity correction. All methods remain research candidates requiring study-specific calibration. Raw PCA reconstruction fractions and FA fitted fractions are expression summaries, not estimates of spatial variance explained.
}
\seealso{\code{\link{passage_test_pathway_gp}}, \code{\link{passage_prepare_gp_kernel}}}
\examples{
set.seed(1)
Y <- matrix(rnorm(60 * 20), 60, 20,
  dimnames = list(paste0("s", 1:60), paste0("g", 1:20)))
xy <- cbind(runif(60), runif(60)); rownames(xy) <- rownames(Y)
fit <- passage_fit_pathway(Y, paste0("g", 1:10), K = 3)
test <- passage_test_pathway_fit(fit, xy, B = 19)
test$p
\dontrun{
gp <- passage_fit_pathway_gp(fit, xy)
gp$summary
}
}
