# PASSAGE

**P**athway **A**ssessment of **S**patial **S**ignatures via **A**ggregated
**G**aussian-process **E**stimation

PASSAGE is an R package for pathway analysis in spatial transcriptomics. It
estimates pathway activity from covariate-adjusted expression and tests for
spatial structure. A separate competitive analysis compares a pathway with
matched gene sets.

Current version: **0.6.1**.

- [Methods](https://diptavo.github.io/PASSAGE/methods/): scientific questions, inference options, intuition and interpretation.
- [Software guide](https://diptavo.github.io/PASSAGE/software/): installation, input format, R functions and command-line usage.
- [Example analysis](https://diptavo.github.io/PASSAGE/example/): analyze three pathways in the prepared kidney sample, including competitive testing.
- [Full prepared cancer datasets](analysis/cancer_data/): 13 complete samples across five tissue groups, all 50 Hallmark pathways, provenance, and download and analysis commands.
- [Interchangeable competitive statistics](docs/competitive_interface.md): global or pathway-local estimation, custom statistics, calibrated combinations, matching diagnostics, and interactive R examples.
- [Release notes](NEWS.md)

## Installation

```r
install.packages("remotes")  # if needed
remotes::install_github("diptavo/PASSAGE")
library(PASSAGE)
packageVersion("PASSAGE")
```

R 4.1 or later and `Matrix` are required. The pathway PCA, factor-analysis and
GP comparison functions need no additional packages. Install `irlba` for fast
global PCA on the full kidney example. `GpGp` and `fields` are optional and are
used by the separate descriptive GP fitting function.

## Prepared kidney example

The repository includes KC2: **2,166 spots, all 17,630 prepared genes**, an
adjustment design, raw counts and three Hallmark pathways: hypoxia,
angiogenesis and inflammatory response. The full gene background supports
competitive testing.

After installing PASSAGE and `irlba`, run this from the repository directory:

```bash
Rscript inst/examples/run_kidney.R
```

It runs pathway PCA, the GP and fixed-kernel self-contained tests, and competitive
testing. Results are written to `results/kidney_all/summary.csv`, with fitted
objects, component scores, loadings and candidate driver genes alongside it.
The default uses three pathway components and 199 null draws for a quick
execution check. To use more draws and a new output directory:

```bash
Rscript inst/examples/run_kidney.R all 9999 results/kidney_9999
```

The [example analysis](https://diptavo.github.io/PASSAGE/example/)
covers installation, running one testing route, PCA versus factor analysis, and
outputs. [Data provenance and licenses](inst/extdata/kidney_KC2_SOURCE.md) describe
the preparation and gene-set snapshot. Execution alone does not establish that
the statistical assumptions hold in this sample.

## Analysis choices

| Question | Estimation | Testing function |
|---|---|---|
| Does spatial covariance improve the model of pathway activity? | PCA or factor analysis within the pathway, followed by component GP fitting | `passage_test_pathway_gp()` |
| Is pathway activity aligned with a prespecified spatial kernel? | PCA or factor analysis within the pathway | `passage_test_pathway_fit()` |
| Does a pathway score more highly than comparable gene sets? | Global factors, pathway-local PCA/FA, or a custom representation | `passage_competitive_test()`; `passage_analyze_bundle()` for multiple pathways |

The first two tests use unsmoothed component scores. The GP comparison fits
spatial variance fractions and ranges and includes their likelihood improvement
in the test statistic. Its null calibration repeats that parameter search.

## Interchangeable competitive statistics

Version 0.6.0 separates statistic calculation from competitive testing. Given
expression `Y`, aligned `coords`, and a vector of pathway gene names `pathway`:

```r
bins <- passage_make_gene_bins(Y, n_mean_bins = 3,
                              n_detect_bins = 2, n_var_bins = 2)
passage_matching_diagnostics(pathway, bins)$supported
engine <- passage_fit_engine_pca(Y, coords, K = 6)
statistic <- passage_make_statistic(
  engine, Y, gene_bins = bins,
  statistics = c("score_z", "score_activity_hotspot", "score_sparse_topk_z")
)
result <- passage_competitive_test(
  pathway, statistic, gene_bins = bins, B = 999, seed = 1,
  combine = "cauchy"
)
result$components
result$combinations
```

The 17 existing competitive statistics and the fitted-correlation statistic
are available through one factory. A custom function of gene indices can use
the same testing procedure. `passage_make_local_gp_statistic()` supplies local
PCA/FA statistics with `sum`, `max`, or `topk` component aggregation, refitting
every matched set. Cauchy, minimum-p, and Fisher combinations are calibrated
against shared matched sets, without independent-p approximations. Prespecify
the statistic panel and primary combiner. The
[complete small-data example](docs/competitive_interface.md) also explains
nuisance covariates, rank sensitivity, Monte Carlo uncertainty, and outputs.

## Quick example

This example creates a small synthetic pathway. It needs no downloaded data.

```r
library(PASSAGE)
set.seed(7)
n <- 100
coords <- cbind(x = runif(n), y = runif(n))
activity <- cbind(sin(5 * coords[, 1]), cos(5 * coords[, 2]))
Y <- activity %*% matrix(rnorm(24), 2, 12) +
     matrix(rnorm(n * 12, sd = 0.5), n, 12)
dimnames(Y) <- list(paste0("spot_", seq_len(n)), paste0("gene_", 1:12))
rownames(coords) <- rownames(Y)

fit <- passage_fit_pathway(Y, colnames(Y), method = "pca", K = 2)
test <- passage_test_pathway_gp(
  fit, coords = coords, kernel_method = "dense", B = 199, seed = 1
)

test$p                 # self-contained pathway p-value
test$component_fits    # spatial fraction, range and variance estimates
fit$scores             # component activity at each spot
fit$loadings           # gene contributions to each component
```

For an analysis of measured expression, provide the adjustment design as `X`
to `passage_fit_pathway()`. Omitting `K` selects the smallest PCA rank explaining
95% of adjusted expression variance. Factor analysis requires an explicit rank:
`method = "factor_analysis", K = 3`.

Use the dense GP calculation for small datasets where feasible. For larger
datasets, `kernel_method = "lowrank"` uses a Nyström approximation and keeps all
spots. Its default 64 features is a starting point; check sensitivity to feature
count. Choose rank rules, kernel settings and the primary test before inspecting
the results.

## Interpretation

PASSAGE provides research methods. The Gaussian rotation and fixed-kernel tests
assume independent Gaussian errors across spots with a common gene covariance,
conditional on a fixed adjustment design. In a small simulation pilot, GP null
rejection rates were 2.0–4.5% at a nominal 5%, but rose to 14.5–25.0% when noise
variance varied spatially. Version 0.6.1 adds a separate, approximate
`passage_test_pathway_gp_wild()` option that preserves spot-specific residual
magnitudes and refits pathway estimation in each bootstrap draw. The original
Gaussian rotation test retains its common-variance assumption.
The new competitive interface also requires study-specific calibration;
interchangeable statistics and empirical combinations do not remove matching
assumptions or establish universal calibration.
For the conditional uniform-gene-label null, use `reference = "universe"`
in competitive testing. Use measured covariates, or prespecified negative
controls through `passage_prepare_nuisance()`, to address nuisance variation.
See the [calibration guide](docs/robust_calibration.md) for assumptions and examples.
The [methods guide](https://diptavo.github.io/PASSAGE/methods/#assumptions)
describes the pilot and its scope.

Single-pathway p-values need correction across the chosen family of pathways
and samples. With `B = 199`, the smallest possible p-value is 0.005; use enough
null draws for the intended testing threshold. GWAS integration is a separate
downstream analysis.

## License

Code: MIT. Copyright 2026 Diptavo Dutta. The bundled KC2 data and Hallmark lists
retain their [CC BY 4.0 source terms](inst/extdata/kidney_KC2_SOURCE.md).
