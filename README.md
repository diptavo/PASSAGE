# PASSAGE

**P**athway **A**ssessment of **S**patial **S**ignatures via **A**ggregated
**G**aussian-process **E**stimation

PASSAGE is an R framework for pathway-level inference in spatial
transcriptomics. It asks whether a predefined gene set has more residual
spatial signal than expected for matched genes, quantifies the strength and
coherence of that signal, and identifies genes that consistently drive it.

## Statistical Scope

PASSAGE separates screening from primary inference.

| Layer | Question | Recommended role |
|---|---|---|
| H1 | Does the pathway contain spatial signal? | Self-contained screen |
| H2 | Does signal remain after measured cell-type adjustment? | Adjusted screen |
| H3 competitive | Is residual pathway signal stronger than in matched random gene sets? | Primary pathway inference |

The competitive H3 null is evaluated by Monte Carlo gene-set resampling after
adjustment for cell-type proportions, technical covariates, and optional
background spatial factors. Matching can account for pathway size, expression,
variance, detection rate, and coexpression module.

The current primary statistic is `score_z`, the standardized spatial
variance-component score. `score_robust_z` is a sensitivity analysis. Raw
`score_Q`, cEPSV, pERSA-derived scores, CSPS, GSPS, HCPS, and other experimental
statistics are retained for benchmarking or descriptive interpretation; they
should not replace `score_z` for confirmatory inference without dataset-specific
null calibration.

## Features

- sparse Vecchia spatial precision construction;
- PCA, spatial-basis, smoothed-PCA, NMF, alternating-GP, and optional TMB
  factor engines;
- self-contained H1/H2 score tests with residual-permutation calibration;
- competitive H3 tests with matched Monte Carlo null gene sets;
- empirical and pathway-size-stratified calibration utilities;
- optional generalized-Pareto tail extrapolation in the large-scale workflow;
- pathway spatial-variance, effective-size, covariance, coherence, hotspot,
  and transferability metrics;
- adaptive top-k driver genes, bootstrap selection frequencies, null driver
  controls, and leave-one-gene-out validation;
- reproducibility workflows for breast cancer, kidney cancer, and DLPFC data.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("diptavo/PASSAGE")
```

Optional engines and calibration paths use packages listed in `Suggests`.
The core implementation requires R 4.1 or newer and `Matrix`.

## Inputs

A standard analysis needs:

1. `Y`: locations by genes normalized expression matrix;
2. `coords`: locations by spatial-coordinate matrix;
3. `pathways`: named list of gene symbols;
4. `Z_CT`: aligned cell-type proportions, when available;
5. additional technical or biological covariates;
6. a background gene universe and matching variables for competitive tests.

Rows must be aligned across expression, coordinates, and covariates. Gene names
must be unique and stored in `colnames(Y)`.

## Quick Start

The bundled simulation runs without external data:

```r
library(PASSAGE)
source(system.file("examples", "run_passage_simulation.R", package = "PASSAGE"))
```

For a self-contained H1/H2 screen:

```r
screen <- passage_run(
  Y = Y,
  coords = coords,
  pathways = pathways,
  Z_CT = cell_type_proportions,
  X = technical_covariates,
  hypotheses = c("H1", "H2"),
  calibration = "permutation",
  n_perm = 999,
  seed = 1
)

head(screen$summary)
```

These p-values answer self-contained questions and can saturate in strongly
structured tissue. Use competitive H3 for primary pathway discovery.

## Competitive H3

Construct the H3 design from cell-type proportions, technical covariates, and
any prespecified background spatial factors. Then compare each pathway with
matched random gene sets:

```r
X_h3 <- cbind(technical_covariates, cell_type_proportions,
              background_spatial_factors)

engine <- passage_fit_engine_pca(
  Y = Y,
  coords = coords,
  X = X_h3,
  K = 6,
  m = 20
)

precomp <- passage_h_precompute(engine, X = X_h3)
gene_bins <- passage_make_gene_bins(Y, gene_names = colnames(Y))

h3 <- passage_conditional_competitive_test(
  engine = engine,
  Y = Y,
  pathway = pathways[["HALLMARK_HYPOXIA"]],
  gene_bins = gene_bins,
  precomp = precomp,
  statistic = "score_z",
  sampler = "module",
  B = 9999,
  seed = 1
)

h3$permutation$p
h3$permutation$enrichment
```

Across pathways, apply FDR to the Monte Carlo p-values. The attainable exact
p-value is `1 / (B + 1)`. Generalized-Pareto extrapolation is intended only for
well-populated, diagnostically acceptable null tails; the empirical p-value
must always be retained.

## Drivers And Interpretation

PASSAGE separates pathway significance from driver selection. After a pathway
passes the competitive test, adaptive top-k selection ranks genes using their
spatial score and pathway contribution:

```r
gene_scores <- passage_gene_score_z(engine, Y, precomp = precomp)

drivers <- passage_sparse_topk_pathway_score_stat(
  engine = engine,
  Y = Y,
  pathway = pathways[["HALLMARK_HYPOXIA"]],
  precomp = precomp,
  gene_scores = gene_scores,
  output = "details"
)

drivers$selected_genes
drivers$score_by_k
```

Single-fit ranks are exploratory. The production driver workflow bootstraps
spots, reports gene selection frequency and rank stability, compares against
matched null pathways, and validates drivers by leave-one-gene-out score loss.

Useful pathway summaries include:

| Metric | Plain interpretation |
|---|---|
| competitive enrichment | Observed pathway score relative to matched null sets |
| `mean_propSV_conditional` | Mean residual spatial fraction across pathway genes |
| ePSV | Scale-weighted residual pathway spatial fraction |
| PC1 spatial fraction | Fraction of pathway spatial covariance concentrated in one coordinated axis |
| mean absolute spatial correlation | Average coherence among fitted spatial gene components |
| effective pathway size | Number of effectively independent genes after correlation |
| transferability | Reproducibility of pathway spatial activity across slices or samples |

cEPSV and related quantities remain useful descriptive metrics even when they
are not calibrated as primary test statistics.

## Data Applications

Reproducible Biowulf workflows are included for breast cancer and kidney
cancer. They download public spatial data, build expression, coordinate,
pathway, covariate, cell-type, and matched-background inputs, submit all
analysis through Slurm, and aggregate competitive H3 and driver results.

```bash
CODE=/data/Dutta_lab/SPATH/PASSAGE
git clone https://github.com/diptavo/PASSAGE.git "$CODE"

# Breast cancer: one Visium sample x two cell-type covariate panels
BREAST_ROOT=/data/Dutta_lab/SPATH/PASSAGE_breast_application
mkdir -p "$BREAST_ROOT"
rsync -a "$CODE/workflows/biowulf/biowulf_cancer_panel/" "$BREAST_ROOT/"
bash "$BREAST_ROOT/submit_pipeline.sh" "$BREAST_ROOT" "$CODE" breast 999 9999

# Kidney cancer: three RCC samples x two cell-type covariate panels
KIDNEY_ROOT=/data/Dutta_lab/SPATH/PASSAGE_kidney_application
mkdir -p "$KIDNEY_ROOT"
rsync -a "$CODE/workflows/biowulf/biowulf_kidney_rcc_gwas/" "$KIDNEY_ROOT/"
bash "$KIDNEY_ROOT/submit_pipeline.sh" "$KIDNEY_ROOT" "$CODE" 999 9999
```

The cell-type inputs currently contain broad marker-derived composition scores,
not probabilistic reference-based fractions. The optional kidney workflow adds
post-hoc MAGMA validation with RCC, clear-cell RCC, and papillary RCC GWAS.
See [docs/data-applications.md](docs/data-applications.md) for prerequisites,
the GWAS command, output locations, and interpretation guidance.

## Repository Layout

```text
R/                         package implementation
inst/tmb/                  optional TMB model
tests/                     simulation and smoke tests
inst/examples/             small runnable simulations
analysis/                  calibration and method-development programs
workflows/biowulf/         dataset and cluster workflows
docs/                      statistical and implementation notes
```

Application instructions are in
[docs/data-applications.md](docs/data-applications.md); statistical calibration
details remain in [docs/calibration.md](docs/calibration.md).

The workflows contain NIH Biowulf SLURM templates. Institutional paths are
defaults from the original analyses, not package requirements. Supply your own
project, reference, GWAS, and output roots when reusing them.

## License

MIT License. Copyright 2026 Diptavo Dutta.
