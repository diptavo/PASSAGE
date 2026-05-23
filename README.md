# PASSAGE

**P**athway **A**ssessment of **S**patial **S**ignatures via **A**ggregated **G**P **E**stimation

PASSAGE is a statistical framework for identifying spatially variable pathways (gene sets) in spatial transcriptomic and proteomic data. It fits a multivariate spatial factor model (Vecchia-approximated LMC), then runs joint variance-components score tests on pre-specified gene sets, producing pathway-level p-values, effect sizes, and a variance-share decomposition.

## Status

**Version 0.1.0 — MVP / development**. All math derived in writing, all functions structurally implemented, but NOT yet validated. See `tests/` for validation scripts that need to be run first.

## What it does

Given spatial expression data $\mathbf{Y}\in\mathbb{R}^{N\times G}$ and a collection of gene sets (e.g. GO BP, max 500 genes each), PASSAGE answers four nested questions per pathway:

| Hypothesis | Question | Adjusts for |
|---|---|---|
| **H1** | Does the pathway have any spatial signal? | nothing |
| **H2** | Does it have spatial signal beyond cell-type composition? | cell-type proportions |
| **H3** | Does it have characteristic spatial pattern beyond the genome-wide background? | cell type + background factor scores |
| **H4** | Does its spatial signature differ across pre-specified regions? | cell type + region indicator |

The headline output is the **variance-share decomposition**:

```
cell-type-explained        =  (Q_H1 - Q_H2) / Q_H1
background-spatial-explained =  (Q_H2 - Q_H3) / Q_H2
pathway-specific spatial   =  Q_H3 / Q_H1
```

Plus three pathway-level effect-size estimators that are equitable across pathways of different sizes:
- `R2_cca` — canonical correlation between observed and fitted spatial signal
- `R2_loo` — leave-one-location-out predictive $R^2$
- `PSVS_range` — spatial-range-weighted proportion of spatial variability

## File layout

```
passage/
├── engine_pca.R                 # Layer 1: PCA two-stage engine (E1)
├── engine_cavi.R                # Layer 1: sparse kernel-ordered CAVI engine (E3)
├── passage_score_h1.R           # Layer 2: H1 (any spatial signal)
├── passage_score_h2.R           # Layer 2: H2 (beyond cell type)
├── passage_score_h3.R           # Layer 2: H3 (beyond background)
├── passage_score_h4.R           # Layer 2: H4 (regional differential)
├── passage_pve.R                # Layer 3: PVE estimators
├── passage_omnibus.R            # Layer 4: cross-H decomposition & reporting
├── passage_summary_stats.R      # Layer 5: summary-statistics fallback
├── run_passage.R                # Top-level driver
├── load_passage.R               # Loader (sources everything in order)
└── tests/
    ├── test_engine_pca.R
    └── test_passage_score_h1.R
```

## Installation

```r
install.packages(c("Matrix", "GpGp", "ucminf", "CompQuadForm"))
# Optional, recommended for large N:
install.packages("sparseinv")  # faster Takahashi selected inversion
```

Then in your R session:

```r
source("load_passage.R")
```

## Usage

### Full joint pipeline

```r
results <- run_passage(
  Y = log_counts,             # N x G expression matrix
  locs = coords,              # N x 2 spatial coordinates
  pathways = gobp_pathways,   # named list of gene-name vectors
  gene_names = rownames(...),
  Z_CT = cell_type_props,     # N x C cell-type proportions, enables H2
  K = 6,
  hypotheses = c("H1", "H2", "H3"),
  adjust_method = "BH"
)

head(results$summary_table)
# pathway | size | p_H1 | p_H2 | p_H3 | cell_type_share | background_share |
# pathway_specific_share | R2_cca | PSVS_range | ...
```

### Summary-statistics fallback (no joint engine fit)

For atlases where you only have per-gene nnSVG outputs:

```r
gene_stats <- nnSVG_to_gene_stats(nnSVG_results)
results <- run_passage_summary(gene_stats, pathways = gobp_pathways)
```

### Single-pathway inspection

```r
engine <- fit_engine_pca(Y, locs, K = 6)
precomp_h1 <- passage_h1_precompute(engine)
res <- passage_test_pathway(engine, Y, my_pathway,
                            Z_CT = cell_type_props,
                            precomp_h1 = precomp_h1,
                            hypotheses = c("H1", "H2", "H3"))
print(res)
# - Per-hypothesis p-values
# - Variance-share decomposition
# - Effect sizes
# - SpASSET-selected spatially-active sub-pathway
```

## Engine choices

| Engine | When to use |
|---|---|
| `fit_engine_pca`  | Fast, deterministic, MVP default. PCA two-stage with per-factor Vecchia GP. |
| `fit_engine_cavi` | Principled. Sparse-loading kernel-ordered LMC with spike-and-slab on A. Slower but yields interpretable sparse loadings. |

Both engines emit the same standardized output triple `(A_hat, theta_hat, D_hat)` plus the Vecchia precision matrices, so the downstream Layer 2/3/4 machinery is engine-agnostic.

## Test variants within each H

Each hypothesis test runs an internal omnibus over:

- **Joint SKAT-style** at multiple weight schemes (`equal`, `var`, `range`)
- **Burden** (single direction $\mathbf{1}_p/\sqrt{p}$)
- **SpASSET subset search** — returns the spatially-active sub-pathway in addition to a p-value

ACAT-Cauchy combination yields one omnibus p-value per H.

## Critical limitations and TODOs

1. **Not yet validated.** Null calibration simulations need to be run before any p-values are trusted. The `tests/test_passage_score_h1.R` script is the starting point.

2. **CAVI engine v_k update** uses a dense matrix inverse instead of true Takahashi selected inversion. Works at Visium scale (N up to ~10k); above that, install `sparseinv` and switch to it.

3. **PVE estimator (c) LOO** currently uses a naive 80/20 train-test split rather than per-location Vecchia conditional predictive. Replace in v2.

4. **H4 calibration** is permutation-based (1000 perms default). For per-pathway speed, an analytical mixed-chi² calibration is on the roadmap.

5. **H3 background fit** is per-pathway-size, cached by size bin. For thousands of pathways this is the main cost; consider precomputing background fits at common sizes.

6. **H2 cell-type adjustment** inherits deconvolution error. Recommended practice: run twice with CARD and cell2location, compare.

## Roadmap

- v0.2: Real Vecchia-LOO predictive R² (replaces 80/20 split).
- v0.2: Takahashi selected inversion in CAVI engine via `sparseinv`.
- v0.3: Analytical H4 calibration; SpatialPCA / NSF / MEFISTO engine wrappers.
- v0.4: Multi-sample extension (slice random effects).
- v0.5: Benchmarking driver vs MEFISTO + projection, nnSVG + ORA, SpatialCorr, dCov-pathway.
- v1.0: Real-data validation paper (DLPFC Visium, Slide-seqV2 cerebellum, breast Visium, METABRIC IMC).

## Mathematical references

- Vecchia, A.V. (1988). Estimation and Model Identification for Continuous Spatial Processes.
- Wu, M.C. et al. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test (SKAT).
- Liu, Y., Chen, S., Li, Z., Morrison, A.C., Boerwinkle, E., Lin, X. (2019). ACAT: A Fast and Powerful p-Value Combination Method for Rare-Variant Analysis in Sequencing Studies.
- Bourotte, M., Heaton, M.J., Banerjee, S. (2024). Computational Considerations for the Linear Model of Coregionalization.
- Bhattacharjee, S. et al. (2012). ASSET: A Subset-Based Approach for Cross-Phenotype Meta-Analysis (the ASSET principle behind SpASSET).
- nnSVG: Weber, L.M., Saha, A., Datta, A., et al. (2023). nnSVG for the scalable identification of spatially variable genes.

## License

MIT.
