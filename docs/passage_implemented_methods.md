# PASSAGE implementation and analysis methods

This document summarizes the method implemented in this codebase during the
10x breast cancer spatial transcriptomics analysis. It is written as a
technical methods record: what model was fit, what hypotheses were tested, what
engines were implemented, what metrics were reported, and how the production
analysis was run on Biowulf.

## 1. Scientific goal

The goal is pathway-level spatial inference for spatial transcriptomics. Instead
of testing one gene at a time and combining many gene-level results, PASSAGE
models the multigene expression matrix for a pathway as a coordinated spatial
signal plus residual expression noise.

The main problem we encountered is that the simple self-contained pathway null,
"no gene in the pathway is spatially variable," is often too easy to reject in
real Visium data. Many genes show at least weak spatial structure because of
tissue architecture, cell composition, library-size gradients, stromal regions,
tumor regions, and inflammatory regions. Therefore, the implementation was
extended beyond self-contained pathway testing to include:

- conditional self-contained tests,
- gene-level SVG-ACAT baselines,
- covariance-aware spatial variance metrics,
- competitive pathway tests against matched background genes,
- pathway coherence tests,
- region-enrichment tests,
- multiple low-rank spatial factor fitting engines.

The working model is:

```text
Y = X B + U A' + E
```

where:

- `Y` is an `n x p` matrix of normalized expression for `n` spots and `p`
  genes in one pathway.
- `X` is an `n x q` design matrix of technical, spatial, cell-composition, or
  background covariates.
- `B` is the `q x p` gene-specific covariate coefficient matrix.
- `U` is an `n x K` matrix of latent spatial factor scores evaluated at spot
  coordinates.
- `A` is a `p x K` loading matrix connecting latent spatial factors to genes.
- `E` is residual gene expression variation not explained by the fitted
  covariates or pathway spatial factors.

The latent factors are intended to capture coordinated spatial pathway
variation. The residual term captures gene-specific and non-spatial expression
variation, including coexpression not explained by the fitted spatial factors.

## 2. Data used

The production analyses were set up for two 10x Genomics Visium breast cancer
datasets:

- `Visium_FFPE_Human_Breast_Cancer`, the 10x Visium FFPE Human Breast Cancer
  DCIS/Invasive dataset.
- `V1_Breast_Cancer_Block_A_Section_1`, the 10x Visium Breast Cancer Block A
  Section 1 dataset.

The single-cell reference for deconvolution was the Wu et al. breast cancer
atlas from GSE176078. The reference was used to estimate spot-level cell-type
proportions, which were then included in the H2 and H3 conditional analyses.

The pathway collections used in production were:

- MSigDB Hallmark gene sets.
- MSigDB KEGG gene sets.

The primary production script is:

```text
scripts/run_10x_msigdb_conditional_passage.R
```

The main Biowulf batch driver is:

```text
jobs_runtime/run_10x_msigdb_all_engines.sbatch
```

The follow-up script for condition-specific competitive testing is:

```text
jobs_runtime/run_10x_msigdb_competitive_condition.sbatch
```

## 3. Preprocessing

The 10x Visium data are read from the filtered feature-barcode matrix and
matched to the tissue-position table. The implemented preprocessing is a
standard Visium-style QC and normalization workflow, not a raw-count analysis.

For each dataset:

1. Read the 10x HDF5 expression matrix and spatial coordinates.
2. Retain in-tissue spots.
3. Compute spot-level QC metrics:
   - total UMI count,
   - number of detected genes,
   - mitochondrial fraction,
   - spatial coordinates.
4. Filter spots using the current defaults:
   - minimum UMI count: 500,
   - minimum detected genes: 200,
   - maximum mitochondrial fraction: 0.25,
   - remove the highest library-size tail above the 0.995 quantile.
5. Collapse features by upper-case gene symbol.
6. Restrict to genes appearing in the selected pathway collection.
7. Retain genes detected in at least `max(20, 1% of spots)` spots.
8. Library-size normalize each spot to 10,000 counts.
9. Apply `log1p` transformation.
10. Scale spatial coordinates to the unit square.

The analysis therefore operates on log-normalized expression, with spot-level
technical and spatial covariates included in the design matrix.

## 4. Covariate designs and nested hypotheses

The implementation uses nested covariate designs to separate different kinds of
spatial pathway signal.

### H1: technical and coordinate-adjusted self-contained spatial signal

H1 tests whether a pathway has residual spatial structure after adjusting for
basic technical and low-order spatial covariates.

The H1 design includes:

- intercept,
- log total UMI,
- log detected genes,
- scaled `x` coordinate,
- scaled `y` coordinate,
- quadratic coordinate terms `x2`, `y2`, and `xy`.

Interpretation:

```text
H1 asks whether the pathway contains spatially structured signal beyond
technical effects and a smooth low-order coordinate trend.
```

This remains a self-contained pathway test. It is useful as a first screen but
is expected to be highly sensitive in real tissue.

### H2: cell-composition conditional spatial signal

H2 adds spot-level cell-type proportions estimated from the Wu et al. breast
cancer single-cell reference.

Interpretation:

```text
H2 asks whether the pathway has spatial signal beyond technical effects,
low-order coordinate trends, and cell-type composition.
```

If H1 is significant but H2 is reduced, the pathway spatial signal may be
largely explained by cell-type composition. If H2 remains strong, the pathway
may vary within cell-type composition strata or may reflect spatial biology not
captured by the deconvolution proportions.

### H3: background spatial factor conditional signal

H3 adds broad background spatial factors, in addition to the H2 covariates.
These background factors are intended to capture global tissue architecture and
large-scale spatial expression programs that are not specific to the pathway.

Interpretation:

```text
H3 asks whether the pathway retains spatial signal after adjusting for
technical effects, cell composition, and broad background spatial patterns.
```

This is stricter than H1 and H2. It is the most useful setting for arguing that
a pathway has pathway-specific spatial structure rather than simply reflecting
global tissue morphology.

### H4: region-enrichment spatial signal

H4 is implemented as a region-enrichment null. The current default creates a
two-level region label from the Wu deconvolution output, using an epithelial
enrichment split. The framework can also accept other two-level annotations,
such as tumor versus stroma, invasive versus DCIS, manual pathology regions, or
image-derived regions.

Interpretation:

```text
H4 asks whether the pathway spatial signal is preferentially enriched across a
defined tissue region contrast.
```

This is not the same as H1 to H3. H1 to H3 ask whether spatial signal remains
after adjustment. H4 asks whether the spatial signal is associated with a
specific region contrast.

## 5. Spatial factor fitting engines

The package now has multiple engines for estimating the latent spatial factor
model. The compatible engines return a common `passage_engine` object, which is
then consumed by score tests, variance metrics, competitive tests, coherence
tests, and region-enrichment tests.

The dispatcher is:

```text
R/passage_fit_engines.R
passage_fit_factor_engine()
```

The comparison helper is:

```text
R/passage_fit_engines.R
passage_compare_factor_engines()
```

### 5.1 Two-stage PCA and Vecchia GP engine

Implemented in:

```text
R/passage_engine_pca.R
passage_fit_engine_pca()
```

Steps:

1. Residualize pathway expression on the selected design matrix `X`.
2. Compute a low-rank PCA representation of the residual matrix.
3. Treat the retained PC score vectors as spatial factors.
4. Fit a spatial covariance model to each factor using a Vecchia approximation.
5. Refit gene loadings and gene-level residual variances.
6. Return factor scores, loadings, residual variances, fitted range parameters,
   and downstream summaries.

This is a two-stage engine. It is not merely "adjusting for PCs." The PCs are
estimated from pathway residual expression and are then modeled as spatial
latent factors. The resulting object contains both the multigene loading
structure and the spatial covariance structure of each factor.

### 5.2 Dynamic rank selection

The PCA engine and compatible engines support fixed-rank and variance-threshold
rank selection.

Current production settings:

```text
rank_method = "variance"
variance_threshold = 0.90
max_K = 12
```

The dynamic rule selects the smallest rank whose cumulative singular-value
variance exceeds 90%, subject to the maximum rank. This reduces sensitivity to a
fixed arbitrary value such as `K = 6`, while still preventing overfitting in
small pathways.

### 5.3 Spatial kernels

The implemented covariance kernels are:

- `matern12`,
- `matern32`,
- `matern52`,
- `exponential`,
- `gaussian`.

The main production runs used:

```text
kernel = "matern32"
```

The Matern 3/2 kernel was used as a default compromise: smoother than an
exponential kernel but less rigid than a Gaussian kernel.

### 5.4 Spatial-basis ridge engine

Implemented in:

```text
R/passage_fit_engines.R
passage_fit_engine_spatial_basis()
```

This engine builds a set of spatial radial-basis and polynomial features from
the spot coordinates, residualizes the basis against `X`, and fits a ridge
smoothed expression surface for the pathway. The smoothed multigene surface is
then compressed into low-rank factors by SVD.

Interpretation:

```text
This engine estimates spatially smooth pathway structure directly from spatial
basis functions, then summarizes that structure as latent factors.
```

It is useful as a fast, stable alternative to PCA-first fitting.

### 5.5 Smoothed PCA engine

Implemented in:

```text
R/passage_fit_engines.R
passage_fit_engine_smoothed_pca()
```

This engine starts from residual PCA factors, then smooths the factor scores
over space using a spatial basis and ridge penalty. Gene loadings are refit
after smoothing.

Interpretation:

```text
This engine keeps the PCA-derived multigene directions but forces the factor
scores to emphasize spatially smooth variation.
```

It often behaves similarly to the spatial-basis engine but remains anchored to
the dominant residual expression axes.

### 5.6 NMF engine

Implemented in:

```text
R/passage_fit_engines.R
passage_fit_engine_nmf()
```

Because log-normalized residual expression can be negative, the residual matrix
is shifted before nonnegative matrix factorization. The NMF scores can then be
spatially smoothed and signed gene loadings are refit.

Interpretation:

```text
This engine searches for parts-based pathway programs instead of orthogonal
PCA directions.
```

In the Hallmark analysis, the NMF engine was usually more conservative than PCA,
spatial-basis, smoothed-PCA, and alternating-GP engines.

### 5.7 Alternating GP engine

Implemented in:

```text
R/passage_fit_engines.R
passage_fit_engine_alternating_gp()
```

This engine alternates between updating factor scores and loadings while
penalizing the factors using spatial covariance structure. It uses a Vecchia GP
approximation to encourage spatially coherent factors.

Interpretation:

```text
This is closer to direct low-rank spatial factor fitting than the two-stage PCA
engine, while remaining much cheaper than a full joint likelihood fit.
```

### 5.8 Joint TMB spatial factor engine

Implemented in:

```text
R/passage_joint_factor_tmb.R
inst/tmb/passage_joint_factor_tmb.cpp
```

The main entry point is:

```text
passage_fit_joint_factor_tmb()
```

This engine fits the model more directly by maximizing a joint likelihood using
TMB with Laplace approximation. It uses:

- latent spatial factor scores,
- Vecchia spatial priors,
- gene-specific residual variances,
- a lower-triangular loading parameterization for identifiability,
- optional nested model comparisons across covariate designs.

Helper functions include:

```text
passage_joint_factor_pve()
passage_select_joint_factor_tmb()
passage_fit_joint_factor_pathway()
passage_fit_joint_factor_pathways()
passage_fit_joint_factor_hypotheses()
```

This engine is conceptually closest to the full model `Y = X B + U A' + E`.
However, it is computationally more expensive. In the production all-engine
Hallmark and KEGG matrix, the compatible engines were run first because the
testing and competitive-metric infrastructure currently consumes the common
`passage_engine` object.

## 6. Self-contained pathway score testing

The original PASSAGE-style score test is implemented around:

```text
R/passage_scores.R
passage_score_test()
```

The workflow wrapper is:

```text
R/passage_run.R
passage_run()
```

The self-contained null is:

```text
H0: the pathway has no spatially structured expression component after
adjusting for the selected covariates X.
```

This is analogous to asking whether the pathway as a set carries any spatial
signal. For real spatial transcriptomics, especially cancer tissue, this null
can be rejected for most pathways. That is biologically plausible but not very
discriminating. This is why the later work focused on conditional tests,
competitive tests, and effect-size metrics.

## 7. Gene-level SVG-ACAT baseline

For each pathway, the pipeline also computes a gene-level spatial variability
baseline:

1. Run an SVG-style test for each gene.
2. Collect gene-level p-values within a pathway.
3. Combine them using the Cauchy combination test, also called ACAT.

The resulting `svg_acat_p` is a gene-level aggregation baseline.

Interpretation:

```text
SVG-ACAT asks whether the pathway contains at least one or more genes with
strong gene-level spatial variability.
```

This is useful as a comparator but has the same practical issue as the
self-contained pathway null: in real tissue, many genes are spatially variable,
so almost all broad pathways can become significant.

## 8. Spatial variance and pathway effect-size metrics

The package implements covariance-aware pathway spatial variance summaries in:

```text
R/passage_competitive.R
passage_pathway_covariance_metrics()
```

These metrics were added because a pathway-level methods paper needs more than
"significant or not significant." The metrics estimate how much coordinated
spatial variation the pathway carries, how concentrated that variation is, and
whether pathway genes move together spatially.

### 8.1 Gene-level mean proportion spatial variance

The simple metric is:

```text
mean_propSV = mean_j spatial_variance_j /
              (spatial_variance_j + residual_variance_j)
```

This is easy to interpret but treats genes as separate units. It averages across
genes and does not fully account for gene-gene covariance or pathway coherence.
Because pathway genes can be coexpressed for reasons unrelated to spatial
structure, `mean_propSV` should be interpreted as a gene-marginal summary, not a
complete pathway covariance metric.

### 8.2 Conditional mean proportion spatial variance

The conditional version is:

```text
mean_propSV_conditional
```

It is computed after residualizing expression on the chosen design matrix. For
H2 and H3, this means it estimates residual spatial variation after adjusting
for cell composition and, for H3, broad background spatial factors.

Interpretation:

```text
mean_propSV_conditional is the average residual spatial fraction across pathway
genes under the chosen conditioning design.
```

### 8.3 Covariance-weighted PVE

The covariance-weighted PVE summaries use spatial and residual covariance
matrices for the genes in the pathway:

```text
S_spatial = covariance of fitted spatial pathway signal
S_residual = shrinkage covariance of residual expression
S_total = S_spatial + S_residual
```

The trace metric is:

```text
cwPVE_trace = tr(S_spatial) / tr(S_total)
```

Generalized PVE values are estimated from the eigenstructure of the spatial
covariance relative to total covariance. The implementation reports summaries
such as:

```text
cwPVE_top
cwPVE_mean
```

Interpretation:

```text
cwPVE_trace is the fraction of total pathway covariance attributable to fitted
spatial structure. cwPVE_top asks whether there is a dominant spatial pathway
axis.
```

These metrics are more pathway-aware than a simple mean across genes because
they operate on covariance matrices.

### 8.4 ePSV: effective pathway spatial variability

The ePSV metric summarizes spatial factor contribution while accounting for the
effective spatial scale of the fitted factor.

Conceptually:

```text
ePSV = weighted spatial factor contribution /
       total pathway variation
```

The spatial factor contribution is based on the fitted loading matrix and
factor variance. Longer-range or more spatially coherent factors receive more
weight than factors that behave like very local noise.

Interpretation:

```text
ePSV estimates how much of the pathway's effective variation is explained by
spatial latent factors, with spatial scale included in the weighting.
```

### 8.5 cEPSV: coherence-weighted effective pathway spatial variability

cEPSV is the more stringent metric introduced in this implementation. It adds a
coherence weight to ePSV.

Conceptually:

```text
cEPSV = sum_k spatial_contribution_k
              x spatial_scale_weight_k
              x pathway_loading_coherence_k
        / total_pathway_variation
```

The coherence term measures whether genes in the pathway load together on the
same spatial factor. A factor that affects only one or two genes strongly but
does not represent a coordinated pathway program receives less weight than a
factor with coherent pathway-wide loadings.

Interpretation:

```text
cEPSV estimates coherent, spatially effective pathway variation. It is not just
the average of gene-level spatial variances.
```

This was added to address the concern that pathway genes can be correlated
beyond spatial correlation. cEPSV favors coordinated pathway-level spatial
programs and downweights isolated gene effects.

### 8.6 Additional diagnostics

The metric function also reports:

- `spatial_eff_rank`, the effective rank of spatial covariance.
- `residual_eff_rank`, the effective rank of residual covariance.
- `residual_mean_cor`, the average residual gene-gene correlation.
- `pc1_spatial_fraction`, the fraction of fitted spatial signal aligned with
  the first pathway expression axis.
- `mean_abs_spatial_cor`, the average absolute correlation in fitted spatial
  signal.

These diagnostics are important for interpreting whether a pathway result is
driven by one dominant spatial axis, many weak spatial axes, or residual
coexpression.

## 9. Competitive pathway testing

Competitive testing was implemented because self-contained tests answer a weak
question for real tissue. A more useful question is whether the pathway has more
spatial structure than comparable background gene sets.

The main implementation is:

```text
R/passage_competitive.R
passage_conditional_competitive_test()
```

The supporting functions are:

```text
passage_gene_spatial_stats()
passage_competitive_gene_stat_z()
passage_competitive_permutation_test()
passage_competitive_fast_context()
passage_fast_cEPSV()
passage_fast_ePSV()
```

### 9.1 Conditional Competitive Null

The Conditional Competitive Null is:

```text
H0: after conditioning on the selected covariates, the pathway's spatial metric
is no larger than expected for matched background gene sets.
```

Matched background sets are sampled from genes with similar expression and
detection properties. This avoids comparing highly expressed, well-detected
pathway genes to poorly measured background genes.

The implemented modes include:

- analytic gene-level Z testing for `mean_propSV_conditional`,
- permutation testing for cEPSV, ePSV, covariance PVE, and coherence-type
  metrics.

The production competitive tests used:

```text
COMP_B = 999
```

### 9.2 Why competitive testing matters

If SVG-ACAT and self-contained PASSAGE are significant for nearly every
pathway, the result is not necessarily wrong. It means the tissue has extensive
spatial biology. However, it is not selective enough for a methods paper.

The competitive null is more defensible:

```text
Instead of asking "is there any spatial signal?", it asks "is this pathway more
spatially organized than comparable gene sets?"
```

This is the inference layer that should be emphasized when ranking pathways and
making biological claims.

## 10. Pathway Coherence Null

The pathway coherence test is implemented in:

```text
R/passage_competitive.R
passage_pathway_coherence_test()
```

The Pathway Coherence Null is:

```text
H0: pathway genes are not more spatially coherent than matched random gene sets.
```

This test focuses on whether genes in the pathway move together in space, rather
than whether individual genes are spatially variable.

Implemented statistics include:

- `pc1_spatial_fraction`,
- `mean_abs_spatial_cor`.

Interpretation:

```text
A significant coherence test supports the claim that the pathway is behaving as
a coordinated spatial program, not merely as a list containing many unrelated
SVGs.
```

## 11. Region-Enrichment Null

The region-enrichment test is implemented in:

```text
R/passage_competitive.R
passage_region_enrichment_test()
```

The Region-Enrichment Null is:

```text
H0: fitted pathway spatial signal is not preferentially enriched across the
specified two-level tissue region contrast.
```

The region label must have two levels. The current default in the 10x breast
cancer analysis is based on an epithelial-enrichment split from the Wu
deconvolution results, but the test is general and can be used with pathology
annotations or image-derived regions.

The statistic is based on region-specific spatial quadratic forms for the
fitted factors. The p-value is obtained by permuting region labels.

Interpretation:

```text
Region enrichment is useful when the scientific question is not just whether a
pathway is spatial, but whether its spatial signal localizes to a biologically
defined tissue compartment.
```

## 12. Inference engines versus fitting engines

The codebase now separates two concepts:

### Fitting engines

Fitting engines estimate the latent pathway spatial model:

```text
Y = X B + U A' + E
```

Implemented fitting engines:

- PCA plus Vecchia GP,
- spatial-basis ridge,
- smoothed PCA,
- NMF,
- alternating GP,
- joint TMB latent factor model.

### Testing engines

Testing engines consume fitted quantities or gene-level p-values:

- self-contained PASSAGE score test,
- SVG-ACAT,
- conditional competitive test,
- pathway coherence test,
- region-enrichment test.

This distinction matters. The same pathway model can be fit by different
engines, and the same fitted object can be evaluated under different null
hypotheses.

## 13. Biowulf production workflow

The all-engine production workflow was configured for:

```text
n_perm = 999
COMP_B = 999
rank_method = "variance"
variance_threshold = 0.90
max_K = 12
K = 6 fallback/default where fixed rank is needed
m = 20 Vecchia neighbors
kernel = "matern32"
engine_n_basis = 60
engine_nmf_iter = 300
engine_alt_iter = 4
```

The engines used in the broad all-engine run were:

```text
pca
spatial_basis
smoothed_pca
nmf
alternating_gp
```

The all-engine production output roots on Biowulf are:

```text
results/passage_10x_all_engines_hallmark_perm999/
results/passage_10x_all_engines_kegg_perm999/
results/passage_10x_all_engines_logs/
```

The production matrix covers:

- two Visium breast cancer datasets,
- Hallmark and KEGG pathway collections,
- H1, H2, H3, and H4 covariate or region-testing settings,
- multiple fitting engines,
- self-contained, SVG-ACAT, competitive, coherence, and region-enrichment
  outputs where applicable.

## 14. Current Hallmark interpretation

The completed Hallmark analyses showed that self-contained tests are highly
sensitive:

- H1 was significant for most Hallmark pathways.
- H2 remained significant for most pathways after cell-composition adjustment.
- H3 remained significant for many pathways after background spatial factor
  adjustment.
- SVG-ACAT was significant for nearly all Hallmark pathways.

This is not surprising. A broad cancer tissue section contains tumor,
epithelial, stromal, immune, vascular, proliferative, and inflammatory spatial
gradients. Many genes and therefore many pathways inherit spatial structure.

The more informative results came from H3 competitive cEPSV testing. Across
engines and both datasets, the robust Hallmark pathways included:

- `HALLMARK_INTERFERON_ALPHA_RESPONSE`,
- `HALLMARK_INTERFERON_GAMMA_RESPONSE`,
- `HALLMARK_ALLOGRAFT_REJECTION`,
- `HALLMARK_COMPLEMENT`,
- `HALLMARK_INFLAMMATORY_RESPONSE`,
- `HALLMARK_IL6_JAK_STAT3_SIGNALING`,
- `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`,
- `HALLMARK_ANGIOGENESIS`,
- `HALLMARK_COAGULATION`,
- `HALLMARK_MYOGENESIS`,
- `HALLMARK_E2F_TARGETS`,
- `HALLMARK_G2M_CHECKPOINT`,
- `HALLMARK_KRAS_SIGNALING_UP`.

Biological interpretation:

- Immune and inflammatory programs show strong spatial organization.
- EMT, angiogenesis, coagulation, and complement likely capture stromal,
  vascular, invasive, and tissue-remodeling structure.
- E2F targets and G2M checkpoint capture proliferative spatial programs.
- KRAS signaling up may capture tumor epithelial or inflammatory-tumor
  interface structure.

Methodological interpretation:

```text
The strongest evidence is not that these pathways are merely significant under
a self-contained null. The stronger evidence is that they remain significant
under H3 competitive cEPSV testing across datasets and multiple fitting engines.
```

## 15. Defending "many pathways are significant"

The implemented analysis should not rely on the claim that all or nearly all
pathways are significant under the self-contained null. That result is expected
in spatial transcriptomics and can invite criticism.

The defensible framing is:

1. Use self-contained PASSAGE and SVG-ACAT as sensitivity screens.
2. State clearly that these screens are expected to be liberal in structured
   tissue.
3. Emphasize conditional inference:
   - H2 asks what remains after cell-composition adjustment.
   - H3 asks what remains after broad background spatial adjustment.
   - H4 asks whether spatial signal localizes to a region contrast.
4. Rank and interpret pathways using effect-size metrics such as cEPSV and
   cwPVE, not only p-values.
5. Use competitive and coherence nulls as the main selectivity layer.
6. Prioritize pathways that replicate across:
   - datasets,
   - fitting engines,
   - related metrics,
   - conditional designs.

This converts the result from "everything is significant" to:

```text
many pathways contain some spatial signal, but a smaller subset shows coherent,
competitive, conditionally robust pathway-level spatial organization.
```

## 16. Limitations

The current implementation is useful but has several limitations that should be
acknowledged.

### 16.1 Pathway gene correlation

Pathway genes can be correlated for non-spatial reasons, such as shared
regulation, cell state, copy number, or technical effects. The simple
`mean_propSV` metric does not fully account for this. The covariance-aware
metrics and cEPSV were added to address this issue, but they are still model-
based summaries and should be validated by simulation and sensitivity analysis.

### 16.2 Background spatial factors

Background factors are intended to capture broad tissue-level spatial patterns.
They may also remove real biological signal if the pathway of interest is part
of a global tissue program. Therefore H3 should be interpreted as a stringent
conditional analysis, not as the only valid analysis.

### 16.3 Deconvolution uncertainty

H2 and H3 treat deconvolved cell-type proportions as covariates. The current
pipeline does not fully propagate uncertainty from the deconvolution step into
the pathway test.

### 16.4 Region labels

The default H4 epithelial-enrichment split is a useful working region
definition, but stronger biological interpretation would require pathology or
image-derived annotations.

### 16.5 Joint likelihood engine

The TMB joint engine is implemented, but the broad all-engine production matrix
currently emphasizes the compatible `passage_engine` engines because the
competitive and coherence tests operate on that common object. A future version
should connect the TMB fit object directly to the same competitive metrics.

## 17. Validation already performed

The following local validation was performed during development:

```text
Rscript tests/test_passage_fit_engines.R
Rscript tests/test_passage_competitive.R
Rscript tests/test_passage_smoke.R
Rscript tests/test_spapath_smoke.R
R CMD check --no-manual --no-build-vignettes
```

The package check passed with warnings and notes related to package hygiene,
optional imports, and bundled result/data paths. The core engine and competitive
tests passed.

Simulation work included:

- null type-I checks,
- 1,000 null replicate runs,
- `n_perm = 999` calibration,
- quick smoke runs on the 10x breast cancer data.

Biowulf validation included:

- a small Hallmark spatial-basis smoke run,
- TMB compilation and wrapper checks under the Biowulf R environment,
- production all-engine Hallmark and KEGG workflows.

## 18. Key files

Core model and engine files:

```text
R/passage_engine_pca.R
R/passage_fit_engines.R
R/passage_joint_factor_tmb.R
inst/tmb/passage_joint_factor_tmb.cpp
```

Testing and metrics:

```text
R/passage_scores.R
R/passage_run.R
R/passage_pve.R
R/passage_competitive.R
```

10x and MSigDB workflows:

```text
scripts/run_10x_msigdb_conditional_passage.R
scripts/deconvolve_10x_with_wu_reference.R
scripts/build_wu_deconvolution_reference.R
scripts/run_10x_svg_acat_pathways.R
scripts/run_passage_competitive_metrics.R
scripts/plot_single_pathway_spatial_heatmap.R
```

Biowulf job scripts:

```text
jobs_runtime/run_10x_msigdb_all_engines.sbatch
jobs_runtime/run_10x_msigdb_competitive_condition.sbatch
```

Tests:

```text
tests/test_passage_fit_engines.R
tests/test_passage_competitive.R
tests/test_passage_joint_factor_tmb.R
tests/test_passage_smoke.R
tests/test_spapath_smoke.R
```

## 19. Manuscript-ready methodological framing

A concise manuscript framing would be:

```text
PASSAGE models pathway expression using a low-rank spatial factor model in
which multigene residual expression is decomposed into coordinated latent
spatial factors and gene-specific residual noise. We fit the model after
conditioning on technical covariates, cell-type composition, and optional
background spatial factors. We evaluate pathway spatial organization using both
self-contained score tests and competitive tests against expression-matched
background gene sets. To avoid relying only on p-values, we estimate
covariance-aware pathway spatial variance metrics, including coherence-weighted
effective pathway spatial variability (cEPSV), which upweights spatial factors
that are long-range and coherently shared across pathway genes. This separates
generic spatial variability from coordinated, pathway-level spatial programs.
```

The most important methodological point is that PASSAGE should be presented as
both a testing framework and an estimation framework. The estimation component
fits coordinated spatial pathway factors. The testing component evaluates those
factors under several null hypotheses, including self-contained, conditional,
competitive, coherence, and region-enrichment nulls.
