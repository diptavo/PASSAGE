\name{passage_test_pathway_gp_wild}
\alias{passage_test_pathway_gp_wild}
\alias{passage_prepare_nuisance}
\title{Nuisance Adjustment and Heteroskedastic Residual Bootstrap}
\description{Construct a measured-covariate design with optional prespecified negative-control factors, or calibrate the pathway GP statistic by refitting on spotwise wild-bootstrap draws.}
\usage{
passage_prepare_nuisance(Y, X = NULL, control_genes = NULL, n_factors = 0L)
passage_test_pathway_gp_wild(Y, pathway, coords = NULL, X = NULL,
  method = c("pca", "factor_analysis"), K = NULL, variance_threshold = 0.95,
  kernel_cache = NULL, ranges = c(0.05, 0.1, 0.2, 0.4),
  fractions = c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97),
  kernel_method = c("lowrank", "dense"), n_features = 64L,
  B = 999L, seed = 1L, aggregation = c("sum", "max", "topk"),
  top_k = 2L, leverage = c("hc2", "none"), keep_null = TRUE)
}
\arguments{
\item{Y}{Finite spots-by-genes expression matrix with unique gene names.}
\item{X}{Adjustment design aligned with expression rows. An intercept is included.}
\item{control_genes}{Prespecified negative-control gene names or indices used only to estimate nuisance factors. Controls must lack the biological signal of interest and should be excluded from the tested gene universe.}
\item{n_factors}{Nonnegative, prespecified number of control PCs to append. Positive values require control_genes.}
\item{pathway}{Pathway gene names or indices.}
\item{coords}{Aligned spatial coordinates; required without kernel_cache.}
\item{method}{PCA or factor analysis, repeated in every bootstrap draw.}
\item{K}{Prespecified pathway rank, or NULL for PCA variance retention. FA requires K.}
\item{variance_threshold}{PCA variance threshold repeated in every bootstrap draw.}
\item{kernel_cache}{Prepared cache with the same coordinates and adjustment design.}
\item{ranges}{Prespecified spatial range grid.}
\item{fractions}{Prespecified spatial variance-fraction grid, including zero.}
\item{kernel_method}{Low-rank or dense spatial kernels.}
\item{n_features}{Maximum number of low-rank kernel landmarks.}
\item{B}{Positive number of bootstrap draws.}
\item{seed}{Nonnegative random seed; the caller's RNG state is restored.}
\item{aggregation}{Sum, maximum or top-k component GP statistic.}
\item{top_k}{Prespecified positive component cap for topk.}
\item{leverage}{HC2 divides each null-model residual row by sqrt(1-h_ii); none omits that correction. Near-unit leverage is rejected.}
\item{keep_null}{Retain null statistics and refitted ranks.}
}
\value{The nuisance helper returns X, control_genes, analysis_genes (controls excluded), leverage/design diagnostics and assumptions. The bootstrap returns the observed GP statistic and component fits, an approximate Monte Carlo p-value, optional null statistics/ranks, configuration and timing.}
\details{
The nuisance helper residualizes control genes against measured covariates, standardizes their residual variances, and appends deterministic left singular vectors. It never uses spatial coordinates or pathway significance to select factors. The number and identity of controls must be chosen before testing. Estimated control factors are not known fixed covariates; contaminated controls can remove true biology. The helper does not infer which genes are valid controls.

The wild bootstrap fits the null mean model X beta separately for pathway genes. It multiplies the entire HC2-adjusted residual gene vector at each spot by one independent Rademacher sign, adds back the fitted mean, then repeats mean adjustment, pathway PCA/FA, rank selection and the complete GP parameter search. Kernel geometry and prespecified choices remain fixed. This retains spot-specific residual magnitudes and joint gene variation while breaking between-spot alignment. It changes calibration, not the observed statistic.

This is an approximate residual bootstrap, not an exact Gaussian rotation test. It assumes independent mean-zero spot errors after a correctly specified mean design. Finite-sample validity is not guaranteed, particularly with high leverage, skewed counts, weakly identified components, spatially dependent errors, or estimated covariates. Negative-control factors are held fixed; their estimation uncertainty is not bootstrapped. Small simulation studies do not establish universal validity. The self-contained null differs from matched-gene-set competitive testing.

The GP profile still uses a constant nugget. Its parameter estimates are descriptive, potentially misspecified summaries under heterogeneous noise; bootstrap calibration does not make them unbiased estimates of a heteroskedastic variance model.
}
\seealso{\code{\link{passage_test_pathway_gp}}, \code{\link{passage_competitive_test}}}
\examples{
set.seed(10)
Y <- matrix(rnorm(30 * 10), 30, 10,
  dimnames = list(paste0("s", 1:30), paste0("g", 1:10)))
xy <- cbind(runif(30), runif(30)); rownames(xy) <- rownames(Y)
adjust <- passage_prepare_nuisance(Y)
ans <- passage_test_pathway_gp_wild(Y, colnames(Y), xy,
  X = adjust$X, K = 2, B = 9, n_features = 10)
ans$p
}
