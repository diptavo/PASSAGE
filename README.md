# PASSAGE

**P**athway **A**ssessment of **S**patial **S**ignatures via **A**ggregated
**G**aussian-process **E**stimation

PASSAGE is an R package for pathway analysis in spatial transcriptomics. It
estimates pathway activity from covariate-adjusted expression and tests for
spatial structure. A separate competitive analysis compares a pathway with
matched gene sets.

- [Methods](docs/method_choices.md): models, hypotheses, calibration and limitations.
- [Software guide](docs/software_guide.md): input format, R examples, command-line usage and outputs.
- [Release notes](NEWS.md)

## Installation

```r
install.packages("remotes")  # if needed
remotes::install_github("diptavo/PASSAGE")
library(PASSAGE)
packageVersion("PASSAGE")   # 0.5.0
```

R 4.1 or later and `Matrix` are required. The pathway PCA, factor-analysis and
GP comparison functions need no additional packages. `GpGp` is optional and
is used only by the separate descriptive GP fitting function.

## Analysis choices

| Question | Estimation | Testing function |
|---|---|---|
| Does spatial covariance improve the model of pathway activity? | PCA or factor analysis within the pathway, followed by component GP fitting | `passage_test_pathway_gp()` |
| Is pathway activity aligned with a prespecified spatial kernel? | PCA or factor analysis within the pathway | `passage_test_pathway_fit()` |
| Is a pathway more spatially structured than comparable gene sets? | Global factor model across the supplied gene universe | `passage_analyze_bundle()` |

The first two tests use unsmoothed component scores. The GP comparison fits
spatial variance fractions and ranges and includes their likelihood improvement
in the test statistic. Its null calibration repeats that parameter search.

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

Version 0.5.0 provides Gaussian-reference research methods. The pathway tests
assume independent Gaussian errors across spots with a common gene covariance,
conditional on a fixed adjustment design. In a small simulation pilot, GP null
rejection rates were 2.0–4.5% at a nominal 5%, but rose to 14.5–25.0% when noise
variance varied spatially. A heteroskedasticity correction is not implemented.
The [methods guide](docs/method_choices.md#validation-and-current-limits)
describes the pilot and its scope.

Single-pathway p-values need correction across the chosen family of pathways
and samples. With `B = 199`, the smallest possible p-value is 0.005; use enough
null draws for the intended testing threshold. GWAS integration is a separate
downstream analysis.

## License

MIT. Copyright 2026 Diptavo Dutta.
