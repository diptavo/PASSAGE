# Competitive testing with interchangeable statistics

PASSAGE 0.6.0 separates **estimation**, **statistic calculation**, and
**competitive calibration**. Estimation can supply factors, loadings, and
spatial kernels. A statistic turns a gene set into one or more numbers.
`passage_competitive_test()` compares those numbers with the corresponding
numbers from matched background gene sets. Custom statistics use the same
testing procedure.

This changes the software organization, not the scientific null: a competitive
test asks whether a pathway scores unusually highly relative to its declared
background. A self-contained test asks whether that pathway has spatial
structure at all. The two questions require different reference distributions.

## A small example in interactive R

This complete example needs no downloads. Its low number of background draws
is for learning the interface; it is not a calibration study.

```r
library(PASSAGE)
set.seed(7)
coords <- as.matrix(expand.grid(x = seq(0, 1, length.out = 8),
                               y = seq(0, 1, length.out = 8)))
Y <- matrix(rnorm(64 * 160), 64, 160)
Y[, 1:12] <- Y[, 1:12] + 0.6 * sin(5 * coords[, 1])
dimnames(Y) <- list(paste0("spot", 1:64), paste0("gene", 1:160))
rownames(coords) <- rownames(Y)
pathway <- paste0("gene", 1:12)

bins <- passage_make_gene_bins(Y, n_mean_bins = 2,
                              n_detect_bins = 1, n_var_bins = 2)
support <- passage_matching_diagnostics(pathway, bins)
support$supported
support$strata

engine <- passage_fit_engine_pca(Y, coords, K = 3, m = 6,
                                range_grid = c(0.1, 0.3), verbose = FALSE)
statistic <- passage_make_statistic(
  engine, Y,
  statistics = c("score_z", "score_activity_hotspot", "score_sparse_topk_z"),
  gene_bins = bins
)
statistic(1:12)                 # numbers only; no p-values

result <- passage_competitive_test(
  pathway, statistic, gene_bins = bins,
  B = 99, seed = 1, combine = "cauchy"
)
result$components              # one competitive p-value per statistic
result$combinations            # calibrated combined competitive p-value
result$matching$strata         # available versus required genes in each bin
result$matching$balance        # observed versus background characteristics
```

The global fit is shared because it does not depend on which set is being
tested. The factory caches residuals and the gene-level quantities needed by
the selected statistics. The default equal-weight `score_z` additionally uses
batched sparse set multiplication. Pathway-specific selections, such as the
best top-gene subset, are repeated for every background set.

## Available global statistics

The factory exposes all 17 legacy competitive choices plus
`mean_abs_spatial_cor`, previously available through the coherence interface.
All are oriented so that larger values mean stronger evidence. Implementation
does not imply that every choice has been validated in every data setting.

| Statistic | Emphasis |
|---|---|
| `score_z` | Standardized loading-weighted spatial quadratic score; default |
| `score_Q` | Unstandardized version of the same quadratic score |
| `score_activity_hotspot` | Spatial structure of the mean standardized pathway expression; directional burden option |
| `score_robust_z` | Winsorized gene scores with a gene-correlation adjustment; sign-insensitive gene evidence |
| `score_sparse_topk_z` | Strongest score over prespecified top-gene subset sizes |
| `score_lse_z` | Enrichment for high gene scores under a two-component mixture |
| `score_factor_hsic` | Association of pathway PCs with the fitted global factors |
| `score_coherence` | Fraction of fitted pathway variation concentrated in its first direction |
| `score_pERSA_z` | Spatial-range and reliability-weighted fitted-variance gene scores |
| `score_pERSA_rank_z` | Rank-standardized version of the preceding gene scores |
| `score_ePSA_rank_z` | Rank-standardized fitted-variance scores without reliability weighting |
| `score_pERSA_active_z` | Enrichment for upper-tail fitted-activity gene scores |
| `mean_propSV_conditional` | Mean fitted-variance fraction across genes |
| `cEPSV` | Spatial-range-weighted fitted-variance fraction with loading coherence |
| `ePSV` | Spatial-range-weighted fitted-variance fraction |
| `cwPVE_trace` | Overall fitted-variance fraction for the pathway |
| `cwPVE_top` | Largest fitted-variance fraction along a pathway expression direction |
| `mean_abs_spatial_cor` | Mean absolute correlation among fitted pathway gene profiles |

Many existing options measure fitted variation or coherence rather than a
direct spatial score. Choose a statistic that matches the intended question.
The word `robust` does not mean a heteroskedasticity correction.

Factory options include `score_weight_scheme`, `covariance`, and `coherence`.
Additional named options are routed to their applicable statistics:
`topk_grid` for the sparse scan; `winsor_q` and `corr` for gene-score aggregation;
`n_pc` for HSIC; `scale_floor_frac`, `active_tail`, and `ridge_frac` for activity
scores; and `shrinkage`/`ridge` for covariance metrics. Unrecognized option names
are errors. For example:

```r
sparse_statistic <- passage_make_statistic(
  engine, Y, "score_sparse_topk_z", topk_grid = c(3, 5, 10)
)
```

## Pathway-local PCA or factor analysis

The local factory fits each pathway separately and calculates GP profile
improvements from its observed, unsmoothed component scores. It repeats that
fit and GP parameter search for **every matched gene set**. Only the
coordinate/design-dependent kernel cache is shared. Competitive calibration
does not require a nested self-contained rotation test for every set.

```r
local_statistic <- passage_make_local_gp_statistic(
  Y, coords, method = "pca", K = 3,
  kernel_method = "dense",
  aggregations = c("sum", "max", "topk"), top_k = 2
)
local_result <- passage_competitive_test(
  pathway, local_statistic, gene_bins = bins,
  B = 99, seed = 1, combine = "minp"
)
local_result$components
local_result$combinations
local_result$estimation[[1]]    # target rank and per-component GP contributions
```

`sum` adds all component improvements; `max` takes the largest; `topk` adds the
largest `top_k` values, using all available components if fewer were retained.
Omit `K` for the PCA variance-retention rule (95% by default). Factor analysis
requires a prespecified rank, for example `method = "factor_analysis", K = 3`.
For larger datasets, use `kernel_method = "lowrank"` and assess sensitivity to
`n_features`. These local competitive options remain research-stage.

The self-contained `passage_test_pathway_gp()` also accepts
`aggregation = "sum"`, `"max"`, or `"topk"` and repeats that aggregation in
every Gaussian reference draw. Choose one aggregation in advance or account
for trying several; its p-values do not answer the competitive question.

## Cauchy, minimum-p, and Fisher combinations

Use `combine = "cauchy"`, `"minp"`, or `"fisher"`. These specify how to combine
component evidence, not an assumed independent theoretical distribution.

For each component, PASSAGE ranks all target/background scores together using
conservative upper-tail ranks. It applies the same combination rule to every
row and then compares the target combined score with background combined
scores. Cauchy endpoints are clipped symmetrically to finite values. There is
no independent-p Fisher chi-square reference and no analytic ACAT p-value in
this competitive interface. Shared background sets retain the dependence
among component statistics during combination.

Choose the primary component panel and combiner before looking at results.
Requesting `combine = c("cauchy", "minp", "fisher")` returns all three for
inspection; it does **not** authorize reporting the smallest resulting p-value
without another selection correction. Empirical combination also cannot repair
a mismatched target/background null.

## Custom statistics and shared background sets

A custom function receives integer gene indices and returns one finite number
or a finite named numeric vector with the same component names on every call.
It must apply the same deterministic calculation and any selection to every
set. The following custom statistic scores a directional map with a fixed
spatial kernel:

```r
Z <- scale(engine$residuals)
spatial_kernel <- exp(-as.matrix(dist(coords))^2 / 0.2^2)
custom_statistic <- function(indices) {
  activity <- rowMeans(Z[, indices, drop = FALSE])
  c(fixed_kernel_burden = as.numeric(crossprod(
    activity, spatial_kernel %*% activity
  )) / max(sum(activity^2), .Machine$double.eps))
}
custom_result <- passage_competitive_test(
  pathway, custom_statistic, gene_bins = bins,
  matched_sets = result$matched_sets
)
```

Supply `alternative = "less"` for a statistic whose smaller values are more
extreme, or a named vector of `greater`/`less` choices for a mixed panel.
Reusing `matched_sets` permits direct comparisons on identical draws. Supplied
sets must contain no duplicates within each set, have the target size, and
preserve its matching strata. The default `reference = "complement"` excludes
target genes. With `reference = "universe"`, target overlap is intentional and
allowed; specify the same reference when reusing sets. Their independence is not
assumed, so Monte Carlo intervals are omitted for externally supplied sets.

## Matching and Monte Carlo diagnostics

Version 0.6.1 adds `reference = "universe"` to the generic test, sampler,
matching diagnostics and bundle workflow. It draws from the complete eligible
universe within strata, including target genes, and tests conditional uniform
gene-set membership. This supplies an exchangeable randomization reference
under that explicit null, even when spot noise or nuisance structure varies.
It does not assert that a curated correlated pathway has uniform labels.
`reference = "complement"` retains the earlier target-excluded comparison and
its calibration limitations. See [the calibration guide](robust_calibration.md).

Run `passage_matching_diagnostics()` before fitting expensive models. Its
`supported` flag and `strata` table show whether every required gene draw is
possible. Optional `gene_bins$module` labels are matched jointly with the
expression bins. Insufficient support is an error, with no silent replacement
sampling or coarsening. Prespecify an appropriate universe and binning scheme;
do not try binning schemes until a preferred p-value appears.

For independently generated sets, each component reports its exceedance count,
number of draws, p-value floor `1 / (B + 1)`, `at_floor` flag, and an exact
binomial interval for the conditional background-tail probability. This
interval describes Monte Carlo uncertainty, not biological uncertainty or
uncertainty in the matching assumptions. It is not constrained to contain the
plus-one p-value estimate. Combined statistics share pooled ranks across all
draws; their exceedance indicators are not independent binomial observations,
so their Monte Carlo interval fields are `NA`.

## Bundle workflow, covariates, and rank sensitivity

`passage_analyze_bundle()` exposes these options for multiple pathways. Its
default remains global fixed-rank estimation and `score_z` with no combination.
Continuing the small example:

```r
bundle <- list(Y = Y, coords = coords,
               pathway_sets = list(example_pathway = pathway))
analysis <- passage_analyze_bundle(
  bundle, B = 99, seed = 1, min_size = 5, gene_bins = bins,
  K = 3, m = 6, range_grid = c(0.1, 0.3), verbose = FALSE,
  statistics = c("score_z", "score_activity_hotspot", "score_sparse_topk_z"),
  combine = "cauchy"
)
analysis$summary               # primary result per pathway
analysis$tests                 # every component and combination
analysis$estimation            # fitted ranks, variance retention, GP settings
analysis$details$example_pathway$competitive
```

Put measured technical covariates in `bundle$X`, with row names aligned to `Y`,
or supply `X` when fitting the global/local estimator directly. The default
global precomputation inherits the engine's design. Explicit precomputation
with a different design remains permitted for legacy compatibility, but should
not be used inadvertently: fit and test with the same nuisance design. Do not
remove a spatial biological effect merely because it explains much variance.

For global automatic rank selection, set `rank_method = "variance"`,
`variance_threshold = 0.90`, and an appropriate `max_K`. Inspect
`rank_info$threshold_reached`: a rank cap may prevent reaching the target.
For a prespecified sensitivity panel, use `rank_grid = c(3, 6, 10)` with
`rank_method = "fixed"`. Each statistic/rank pair becomes a separate component;
a requested combined test includes them all. Without combination, the first
component remains the primary result. More factors are not guaranteed to
remove nuisance sensitivity or improve power.

For bundle pathway-local estimation, use `estimation = "pathway"`,
`local_K = 3`, and `gp_aggregations = c("sum", "max")`; keep `statistics` at its
default because `gp_aggregations` selects the local statistics. A custom
function can instead be supplied as `statistics`; it must index the same
nonconstant gene universe retained by the bundle.

The bundle selects the first requested combination as primary, or the first
statistic when no combination was requested. `summary$q_BH`/`q_BY` adjust that
primary family across pathways. In `tests`, adjustment is performed across
pathways separately for each named test. Neither adjusts selection among
statistics or combiners after seeing results.

## Scope of the changes

Scoring now applies the global engine's gene centers and scales consistently
with its fitted loadings and residual variances. Pathway-holdout engines now
reset preprocessing metadata to the unscaled units used in their final
all-gene loading refit. Previously scaled analyses should be rerun.

The new interface, statistic choices, and combination calibration do not
establish universal type-I-error control. Competitive matching must suit the
study, and local PCA95 showed an iid-null calibration warning in the earlier pilot.
Spatial heteroskedasticity remains a known limitation of the self-contained
Gaussian reference. Retain study-specific null simulations and independent
validation before relying on a new statistic/default in substantive analyses.

For large pathway collections, `passage_analyze_bundle()` defaults to
`keep_matched_sets = FALSE`: it retains scores, matching summaries and seeds,
but avoids storing every sampled gene index for every pathway. Set it to TRUE
when explicit set membership is needed. The lower-level
`passage_competitive_test()` retains sets by default. Local diagnostics retain
rank and variance retention for all sets, and detailed component fits for the
observed pathway.

## Validation of version 0.6.0

The [October 30 validation report](validation/competitive_20261030/findings.md)
contains the completed 1,200-dataset pilot, uncertainty intervals, and software
checks. Directional activity improved dense-signal detection, but combined
null rejection rates were approximately 5–8% at a nominal 5% and individual
statistics showed calibration warnings under heteroskedasticity or nuisance
variation. The default remains unchanged; the combinations are optional.
The main pilot covers global statistics. Local GP aggregation passed smoke
and regression checks, but its power and calibration need separate evaluation.
