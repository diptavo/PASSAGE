# Unequal noise and nuisance variation

Version 0.6.1 adds three explicit options. They retain the estimation/statistic/
testing separation. Choose them before looking at pathway significance.

## Competitive testing: full-universe randomization

```r
engine <- passage_fit_engine_pca(Y, coords, X = X, K = 4, verbose = FALSE)
bins <- passage_make_gene_bins(Y, n_mean_bins = 3,
                              n_detect_bins = 2, n_var_bins = 2)
statistic <- passage_make_statistic(
  engine, Y, statistics = c("score_z", "score_activity_hotspot"),
  gene_bins = bins
)
result <- passage_competitive_test(
  pathway, statistic, bins, B = 999, seed = 1,
  reference = "universe", combine = "minp"
)
result$components
result$combinations
result$matching$mean_target_overlap
```

Here `Y` has spots in rows and genes in columns, `coords` has the same spot
order, `X` is a numeric covariate matrix (or NULL), and `pathway` contains gene
names. The bundle API accepts the same `reference` argument.

The older default, `reference = "complement"`, excludes the target from every
background set. This is a target-versus-complement comparison and lacks the
following exact conditional randomization argument in a finite gene universe.

For `"universe"`, condition on expression, covariates, matching strata and any
shared global fit. Suppose the target is uniform among all gene sets with its
size and stratum counts. Draw each reference set independently from that same
space and apply the identical statistic. The observed and reference scores
are then exchangeable. The plus-one Monte Carlo tail is super-uniform; the
symmetric pooled-rank combinations preserve this property. This argument
does not assume homoskedastic spots or a correctly specified global GP model,
because expression is held fixed under this **gene-label null**.

Any shared fit, design or stratum definition must treat target designation
symmetrically. A global PCA fitted once without reference to the target meets
this computational requirement. A model fitted specifically by holding out
the observed target does not automatically meet it: repeat that set-dependent
fitting rule for every reference set instead. The statistic callback must not
privilege the originally named pathway.

Target overlap is part of the null. Removing overlapping draws destroys this
construction. Genes never repeat within a set, but sets can recur across
draws. If strata leave only one possible set, p-values are one and diagnostics
report no reference variation. Many small strata can make the test uninformative.

This is not a proof for a fixed curated pathway under repeated expression
sampling. Co-regulated modules can violate uniform membership within bins;
module-aware nulls and independent calibration remain necessary. More draws
improve Monte Carlo precision but cannot repair a wrong scientific null.

## Nuisance variation: use a declared mean design

Measured technical variables are preferred:

```r
X <- model.matrix(~ batch + log_library_size, data = spot_metadata)
engine <- passage_fit_engine_pca(Y, coords, X = X, K = 4, verbose = FALSE)
```

Use covariates relevant to the scientific question. Cell-type proportions,
for example, change the interpretation to spatial variation beyond composition.
Coordinates alone should not be regressed out indiscriminately: that can remove
the spatial biology being tested. A larger global K does not by itself remove
nuisance; it can allocate more factors to it.

When validated negative-control genes exist, their PCs can augment measured
covariates without using pathway significance or target-gene expression:

```r
adjustment <- passage_prepare_nuisance(
  Y, X = X, control_genes = negative_controls, n_factors = 3
)
X <- adjustment$X
Y <- Y[, adjustment$analysis_genes, drop = FALSE]
adjustment$diagnostics
# Rebuild bins, engine and statistic using this Y and X.
```

Controls must share unwanted variation but lack the biology of interest.
Neither this helper nor PCA can certify that assumption. Prespecify controls
and rank; exclude controls from the tested gene universe. This can lose real
signal if controls are contaminated. The returned design includes leverage
and residual degrees of freedom. Control-factor estimation uncertainty is not
accounted for by treating those columns as fixed in self-contained inference.

## Self-contained GP testing: preserve unequal noise

```r
result <- passage_test_pathway_gp_wild(
  Y, pathway, coords, X = X, K = 2,
  B = 999, seed = 1, aggregation = "sum"
)
result$p
result$component_fits
result$max_leverage
```

The observed PCA/FA and GP statistic are unchanged. Each null draw multiplies
the entire residual gene vector at a spot by one random sign, after HC2 leverage
adjustment. It then repeats covariate adjustment, pathway estimation, rank
selection and the GP range/fraction search. The same sign across genes preserves
their joint residual variation; keeping magnitudes at their original spots
retains the noise gradient that Gaussian rotations erase.

This is an **approximate residual wild bootstrap**, not an exact test. It
requires a suitable mean model and independent spot-error vectors. It does not
solve spatially correlated residual errors, omitted nonlinear means, contaminated
controls or arbitrary count-distribution effects. It holds the supplied design
fixed. High leverage and weakly separated PCs warrant particular care.

The original `passage_test_pathway_gp()` remains the Gaussian rotation option.
The wild method is a self-contained test; it does not replace competitive
backgrounds. Parameter estimates remain descriptive plug-in fits with no
bootstrap confidence intervals. The fitted GP still uses a constant nugget:
under heteroskedasticity, its fitted spatial fraction/range should not be read
as unbiased physical variance parameters. This option adjusts the null
calibration; it does not fit a heteroskedastic GP variance surface.

## Validation

The [completed October 31 report](validation/calibration_20261031/findings.md)
compares both competitive references on the same datasets, with measured and
idealized negative-control nuisance designs, and compares Gaussian rotations
with wild-refit calibration on the same pathway statistics. The main experiment
contains 6,800 datasets; a separately labeled weak-signal follow-up contains 600.
All remote numerical work used Biowulf batch jobs with independent scenario
seeds and 399 draws. Both audits found no computational failures or warnings.

Full-universe minimum-p rejection was 4.1%, 4.9% and 4.8% under unequal noise,
strong nuisance and both, respectively. This combination used the loading,
directional, robust and sparse panel. Wild-bootstrap rejection was 4.2% under
unequal noise versus 14.4% with rotations, and 6.2% with noise plus nuisance
versus 12.4%. The combined case remains uncertain (wild 95% CI 4.3–8.7%).
Three adjusted robust/sparse competitive checks had pointwise calibration
warnings at 7.2–7.4%; these are retained in the report. Neither API is advertised
as universally calibrated. The follow-up showed essentially unchanged weak
signal detection with common variance; lower detection under unequal variance
must be interpreted alongside the correction of inflated null rejection.

Statistical background: Goeman and Buehlmann, Bioinformatics 2007,
doi:10.1093/bioinformatics/btm051; Davidson and Flachaire, Journal of Econometrics
2008, doi:10.1016/j.jeconom.2008.08.003. These motivate the null distinction and
wild-bootstrap construction; they do not establish a theorem for PASSAGE's
adaptive GP statistic.
