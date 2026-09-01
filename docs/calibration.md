# Calibration And Inferential Recommendations

## Primary Null

The primary PASSAGE claim is competitive H3:

> After adjustment for measured cell-type composition, technical covariates,
> and prespecified background spatial structure, does this pathway have more
> residual spatial signal than a matched random gene set?

This is stronger and more selective than the self-contained H1/H2 question of
whether any residual spatial signal exists in the pathway.

## Recommended P-Value Construction

1. Fit or define the H3 adjustment without using the tested pathway to select
   covariates.
2. Compute the pathway `score_z` statistic.
3. Sample random gene sets matched on pathway size and gene-level nuisance
   features. Include coexpression-module matching when feasible.
4. Recompute the complete statistic for every random set.
5. Use `(1 + number(null >= observed)) / (B + 1)`.
6. Apply BH FDR to the resulting competitive p-values across pathways.
7. Repeat with `score_robust_z` and alternative reasonable matching designs as
   sensitivity analyses.

For extremely small p-values, generalized-Pareto extrapolation may supplement
the empirical tail only when the threshold, exceedance count, shape, scale,
and fit stability are acceptable. Report the empirical Monte Carlo bound and
whether extrapolation was used.

## Current Empirical Reading

- `score_z` has shown the most reliable null behavior across the breast,
  kidney, and DLPFC calibration settings tested to date.
- `score_robust_z` is useful when a few genes dominate, but remains a secondary
  sensitivity statistic.
- Raw `score_Q` is more dependent on pathway size and covariance structure;
  standardization or empirical null calibration is needed.
- cEPSV is interpretable as a descriptive effect-size concept, but its null
  behavior was not adequate as the primary resampled test statistic.
- pERSA-derived and newer CSPS/GSPS/HCPS statistics remain experimental until
  larger null and power benchmarks establish stable operating characteristics.
- Self-contained H1/H2 can be highly saturated in spatially heterogeneous
  tissue. This is expected behavior for their null and is not fixed by merely
  changing the analytical tail approximation.

## Required Study-Specific Checks

Before confirmatory analysis, stratify null calibration by pathway size and,
where relevant, expression, detection rate, gene correlation, spatial range,
and cell-type-reference choice. Report type-I error with uncertainty at the
target alpha levels, Monte Carlo resolution, sampler acceptance or effective
sample size, and the number of pathways in each stratum.

Driver-gene claims require bootstrap selection frequencies, matched-null driver
frequencies, and leave-one-gene-out or held-out-slice validation. A single
adaptive top-k list is not a calibrated driver result.
