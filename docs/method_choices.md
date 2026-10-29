# Methods

PASSAGE separates estimation of pathway activity from the question being
tested. A self-contained test asks whether the retained pathway components
have spatial structure after adjustment. A competitive test asks whether the
pathway has more spatial signal than comparable genes. These are different
null hypotheses; significance under one does not imply significance under the
other.

## Expression and adjustment

Let $Y$ be a spots-by-genes normalized expression matrix and $X$ a full-rank
adjustment design containing an intercept. For a pathway $P$, estimation starts
from $Y_P$ after projection off $X$.

The covariates define the interpretation. Adjusting for technical measurements
asks about signal beyond those measurements. Adding cell-composition estimates
asks about signal beyond the estimated composition. Composition adjustment can
remove biological variation of interest, so it should follow the scientific
question rather than serve as an automatic default.

## Pathway estimation

`passage_fit_pathway()` uses only the selected pathway genes. It does not use
coordinates or a previously fitted global model. Missing genes and genes with
zero adjusted variance are recorded in the result.

| Estimator | Definition | Rank selection |
|---|---|---|
| PCA | SVD of adjusted expression, without scaling each gene to unit variance | Smallest rank explaining 95% of variance, or a prespecified `K` |
| Factor analysis | Maximum-likelihood FA of the gene correlation matrix, with unrotated regression scores | An explicit, identifiable `K` |

Both estimators return spot scores and gene loadings. FA requires a nonsingular
gene correlation matrix; it is unsuitable when the available residual degrees
of freedom cannot support that matrix and the requested factor count. All
subsequent tests concern the retained representation. A weak spatial component
discarded during estimation may be missed.

## GP covariance comparison

`passage_test_pathway_gp()` fits a spatial covariance model to each unsmoothed
component. In an orthonormal covariate-residual basis, its working covariance is

$$
C_k = \sigma_k^2\{(1-h_k)I + h_k K_{\ell_k}\},
$$

where $K_{\ell_k}$ is the projected Matérn 3/2 kernel at range $\ell_k$, normalized
to have trace equal to the residual degrees of freedom. Thus $h_k$ is the
average spatial variance fraction in that residual space. The null is $h_k=0$
for every retained component.

The implementation profiles out $\sigma_k^2$ and maximizes the profile
restricted log likelihood over prespecified range and fraction grids. The
pathway statistic sums twice the component improvements over their nonspatial
models. The fitted ranges and fractions therefore contribute directly to
testing. For correlated FA scores, this sum is a composite statistic, not a
joint independent-factor likelihood ratio.

Calibration uses joint Gaussian rotations of the scores in the residual space.
Each draw preserves the complete score Gram matrix and repeats the same
range/fraction search. A chi-square reference is not used. The Monte Carlo
p-value is $(1 + \#\{T_b\geq T_{obs}\})/(B+1)$.

Under the fixed-design Gaussian null, pathway PCA and the implemented FA
estimator depend on the gene Gram matrix and transform with row rotations.
Their loadings, rank rule and score Gram matrix can therefore be conditioned
on. There is no need to refit PCA or FA numerically for every rotation; the GP
parameter search is repeated because it changes under rotation.

### Kernel calculation

The default approximation uses coordinate-selected Nyström landmarks. A shared
basis and cached kernel spectra allow batched likelihood calculations without
forming full-size null score matrices. `passage_prepare_gp_kernel()` can reuse
this work across pathways with exactly the same spots, order and design.

The dense option provides a reference for smaller datasets, with a default
limit of 1,000 spots. Low-rank kernels retain all spots but approximate their
covariance. Feature count can affect the statistic and power; 64 features is
not an accuracy guarantee. Compare prespecified feature counts and use a dense
reference where affordable. Do not choose the approximation giving the smallest
p-value.

Reported parameters are grid estimates without confidence intervals. A zero
spatial fraction leaves the range unidentified. Boundary flags indicate that
an estimate lies at an edge of the specified grid.

## Fixed-kernel subspace tests

`passage_test_pathway_fit()` orthonormalizes the adjusted score span to obtain
$U$. For a fixed, positive-semidefinite spatial kernel $R$ and residual degrees
of freedom $d$, the two statistics are

$$
T_{trace}=d\,\mathrm{tr}(U^\top R U),\qquad
T_{root}=\lambda_{max}(U^\top R U).
$$

The trace summarizes alignment across the retained span; the largest root
emphasizes its strongest spatial direction. Both are invariant to invertible
changes of component coordinates. Their null reference simulates Gaussian
subspaces in the covariate-residual space. An optional within-pathway omnibus
p-value applies Bonferroni correction to the two tests.

The default kernel is Matérn 3/2 with a prespecified range. A range estimated
from the same scores cannot simply replace it without calibrating that
selection. This dense implementation has a default 3,000-spot limit and never
silently subsets the data.

The scale-free largest-root test is not an ordinary Tracy–Widom test. PC-wise
SKAT with Davies calibration and its Fisher or ACAT combinations are not
implemented by these pathway APIs.

## Competitive analysis

`passage_analyze_bundle()` fits a global PCA factor engine across the supplied
gene universe and uses its loadings, residual variances and spatial kernels to
score each pathway. It compares the observed `score_z` with scores of matched
gene sets. The p-value is an empirical upper-tail probability, not a
standard-normal tail despite the statistic's name.

Default matching uses expression, detection and variance bins. An optional
module column adds module matching. Each draw preserves stratum counts, contains
no duplicate genes and excludes the tested genes. Inadequate background support
causes an error rather than a fallback to unmatched genes. Sets may recur
across separate Monte Carlo draws.

This procedure relies on exchangeability within the chosen matching strata;
matching a few gene summaries does not establish that assumption for every
dataset. The API returns BH and BY adjustments across the eligible pathways in
one call. Multiple samples may require a larger prespecified testing family.
Driver rankings are exploratory and need stability and independent validation.

This global procedure is not a competitive version of the pathway GP test.
A pathway-local competitive GP analysis would need to repeat pathway estimation
and GP fitting for each matched set and is not currently implemented.

## Descriptive fits and downstream analysis

`passage_fit_pathway_gp()` optionally fits each component using GpGp and reports
covariance estimates and convergence. It is a separate descriptive API. Its
outputs are not passed into either pathway test, and it does not report spatial
test p-values. The GP comparison fits its own covariance parameters internally.

GP-smoothed predictions inherit spatial structure from the smoother, even under
a nonspatial expression null. The testing functions therefore use unsmoothed
estimated scores.

GWAS integration follows the spatial analysis. It requires a separately defined
gene or pathway hypothesis, appropriate variant-to-gene mapping and an LD-aware
association or enrichment analysis. A PASSAGE run does not perform those steps
or supply genetic evidence by itself.

## Validation and current limits

The pathway references assume independent Gaussian row errors with a common
gene covariance, conditional on a fixed adjustment design. Arbitrary
correlation between genes is allowed. Normalized counts need not meet this
model. Library-size or composition covariates estimated from the same expression
are not automatically externally fixed for the calibration argument.

A pilot used 100 spots on a 10 × 10 grid, 12 genes, an intercept and one fixed
nuisance covariate, with 200 replicates per scenario. Correlated-gene scenarios
used three orthonormal loading directions with variances 4, 3 and 2, and gene
noise variance 0.08. Alternatives assigned a spatial fraction of 0.2 to one or
all three factors, using Matérn 3/2 range 0.2. The GP test used 32 landmarks,
four candidate ranges and 199 null draws. Rejection rates at a nominal 5% were:

| Estimation | Independent-gene null | Correlated-gene null | One spatial factor | Three spatial factors | Heteroskedastic null stress |
|---|---:|---:|---:|---:|---:|
| PCA, 95% variance | 4.0% | 2.0% | 21.5% | 69.0% | 25.0% |
| PCA, K = 2 | 2.5% | 3.5% | 28.0% | 73.0% | 14.5% |
| FA, K = 2 | 4.5% | 3.5% | 27.5% | 74.5% | 14.5% |

All six Gaussian-null binomial 95% intervals included 5%. With 200 replicates,
the pilot gives only a preliminary check, especially at smaller significance
thresholds. In the heteroskedastic scenario, each independent spot's entire
random expression vector was multiplied by a scale proportional to
$\exp\{1.5(x_i-0.5)\}$, normalized to have mean squared scale one. This creates
about a 20-fold variance gradient without cross-spot correlation and violates
the common-variance null. The resulting inflation is a material robustness
failure. The fixed-kernel tests also inflated under that stress.

A variance-aware null and its calibration remain development work. Current
p-values should not be treated as broadly calibrated for heterogeneous spatial
transcriptomics data. Choose the estimator, kernel settings and testing family
in advance; assess calibration and approximation sensitivity; then interpret
replicated results with the covariate adjustment made explicit.
