# Introduction

Spatial transcriptomics measures gene expression while preserving the physical
organization of tissue. This spatial context is essential for studying cancer,
where malignant, stromal, immune, vascular, and epithelial compartments are
organized into structured neighborhoods rather than randomly mixed cell
populations. In breast tumors, spatial localization of transcriptional programs
can reflect tumor invasion, immune infiltration, stromal remodeling, epithelial
differentiation, angiogenesis, and proliferative architecture. Methods that
ignore tissue coordinates can identify differentially expressed genes or
pathways, but they do not directly ask whether a biological program is organized
in space.

Most spatial transcriptomics analyses begin with gene-level spatially variable
gene discovery. These approaches test whether each gene exhibits spatial
autocorrelation, spatial smoothness, or excess spatial variance across tissue
locations. Gene-level spatial testing is useful for feature discovery and
visualization, but it has limitations for pathway-level biology. First, pathway
activity is typically distributed across multiple genes, each of which may have
modest spatial signal. Second, gene-level tests can be sensitive to technical
variation, library-size gradients, cell-type composition, and broad tissue
architecture. Third, aggregating individual spatially variable genes does not
distinguish a coordinated pathway program from a gene set that simply contains
many unrelated spatially variable genes. These issues are especially important
in cancer tissues, where many genes can be spatially variable because of strong
histologic structure.

Gene set and pathway analyses provide a natural way to move from individual
genes to interpretable biological programs. However, standard gene set analyses
were not designed for spatial transcriptomics. Many pathway scoring strategies
summarize expression across genes and then compare scores between groups, while
many enrichment strategies test whether a list of genes is overrepresented in a
predefined pathway collection. These approaches can identify pathways associated
with a phenotype or region, but they do not directly model the spatial
covariance structure of a pathway across tissue. A pathway-level spatial method
should account for the fact that genes within a pathway may be correlated, that
pathway signal may be low-rank, and that spatial structure may remain after
adjustment for technical covariates, cell-type composition, and broad background
spatial patterns.

We developed PASSAGE as a pathway-level spatial inference framework for spatial
transcriptomics. PASSAGE models the multigene expression matrix for a pathway as
the sum of observed covariate effects, coordinated latent spatial factors, and
residual expression variation. This formulation treats a pathway as a
multivariate spatial object rather than as a collection of independent genes.
The low-rank structure provides computational efficiency and reflects the
biological expectation that many pathway genes may share a smaller number of
spatial expression programs. The spatial factor formulation also separates
estimation from inference: first, a pathway-level spatial representation is
estimated; second, the fitted representation is evaluated under several null
hypotheses.

The first inferential target is a self-contained spatial pathway null: whether
a pathway contains any residual spatial signal after adjustment for a selected
covariate design. Although this test is useful as a sensitive screen, it can be
too permissive in structured tissues, because many broad pathways contain at
least one spatially variable gene. PASSAGE therefore implements additional
inferential layers. Conditional tests evaluate pathway spatial signal after
adjusting for cell-type composition and background spatial factors. Competitive
tests compare pathway spatial metrics against expression-matched background
gene sets. Coherence tests ask whether pathway genes show greater coordinated
spatial structure than matched gene sets. Region-enrichment tests ask whether
the fitted spatial pathway signal is enriched across a specified tissue-region
contrast. Together, these tests allow the analysis to move beyond the question
of whether any pathway gene is spatially variable and toward the more specific
question of whether a pathway forms a coherent, conditionally robust spatial
program.

PASSAGE also estimates pathway-level spatial effect sizes. In addition to
gene-marginal summaries such as the mean proportion of spatial variance, we
implemented covariance-aware pathway spatial variance metrics. These include
trace-based covariance-weighted spatial variance, generalized pathway spatial
variance summaries, effective pathway spatial variability, and a
coherence-weighted effective pathway spatial variability metric. The
coherence-weighted metric is designed to emphasize spatial factors that are both
spatially effective and shared across pathway genes, thereby reducing reliance
on isolated gene-level effects.

Here, we describe the PASSAGE model, fitting engines, preprocessing workflow,
covariate designs, pathway hypotheses, spatial variance metrics, and inference
procedures implemented for breast cancer Visium analyses using MSigDB Hallmark
and KEGG pathways. The methods are implemented in R and were designed to run
efficiently on high-performance computing infrastructure while retaining
multiple sensitivity engines for fitting the latent spatial factor model.

# Methods

## Overview

PASSAGE is a pathway-level spatial modeling framework for spatial
transcriptomics. For each pathway, PASSAGE fits a low-rank spatial factor model
to normalized expression across spatial spots. The fitted model is then used to
perform self-contained spatial pathway tests, gene-level spatially variable gene
aggregation, competitive pathway tests, pathway coherence tests,
region-enrichment tests, and pathway-level spatial variance estimation.

The implemented analysis workflow consists of the following stages:

1. Spatial transcriptomics preprocessing and quality control.
2. Pathway gene-set loading and gene matching.
3. Single-cell reference processing and spot-level cell-type deconvolution.
4. Construction of technical, spatial, cell-type, and background-factor
   covariate designs.
5. Fitting of low-rank spatial pathway factor models.
6. Self-contained pathway testing under nested covariate designs.
7. Gene-level spatial testing and ACAT pathway aggregation.
8. Estimation of pathway spatial variance and coherence metrics.
9. Competitive and coherence testing using matched background gene sets.
10. Region-enrichment testing for specified two-level tissue contrasts.

All analyses were implemented in R. Computationally intensive breast cancer
analyses were configured for Biowulf high-performance computing jobs.

## Spatial transcriptomics data

The workflow was prepared for two publicly available 10x Genomics Visium breast
cancer datasets:

- 10x Visium FFPE Human Breast Cancer DCIS/Invasive.
- 10x Visium Breast Cancer Block A Section 1.

For each dataset, the filtered feature-barcode HDF5 matrix and Visium spatial
position files were used as input. Barcodes were matched between the expression
matrix and the spatial tissue-position file, and only spots annotated as
in-tissue were retained before quality control.

## Spatial transcriptomics preprocessing

Raw counts were read from the 10x HDF5 matrix. Gene names were converted to
upper-case symbols, and duplicated symbols were collapsed by summing counts
across features mapping to the same symbol. Spots were filtered using library
size, detected gene count, mitochondrial fraction, and upper-tail library-size
criteria.

The implemented spot-level quality-control thresholds were:

- total UMI count at least 500,
- detected genes at least 200,
- mitochondrial fraction at most 0.25,
- total UMI count below the 0.995 empirical quantile.

Mitochondrial genes were identified using the `MT-` prefix. After spot-level
filtering, genes were restricted to those appearing in the selected pathway
collection. A gene was retained if it was detected in at least
`max(20, 0.01 n)` spots, where `n` is the number of retained tissue spots.
Genes with near-zero variance after normalization were removed.

Counts were library-size normalized to 10,000 counts per spot and transformed
using `log1p`. Spatial pixel coordinates were scaled to the unit square for
model fitting. The resulting expression matrix is denoted by
`Y`, with rows corresponding to spatial spots and columns corresponding to
genes.

## Pathway definitions

Pathways were loaded from MSigDB gene-set files in GMT format. The analyses were
configured for:

- Hallmark gene sets,
- KEGG gene sets.

Pathway genes were converted to upper-case symbols and intersected with genes
available after spatial transcriptomics preprocessing. Pathways were retained if
their post-filtering size fell within the configured pathway-size range. The
default minimum pathway size was 5 genes and the default maximum pathway size
was 500 genes.

## Single-cell reference and cell-type deconvolution

The Wu et al. breast cancer single-cell atlas from GSE176078 was used as a
reference for cell-type deconvolution. A coarse pseudobulk reference matrix was
constructed from the single-cell atlas, with rows corresponding to genes and
columns corresponding to coarse cell types. Reference marker genes were selected
for each cell type using a marker score defined as the expression of the gene in
that cell type minus the maximum expression of the same gene across all other
cell types. Marker genes were retained if they exceeded a minimum marker score
and were expressed in the corresponding reference cell type. The default
configuration selected up to 150 marker genes per cell type.

Spot-level cell-type proportions were estimated using a constrained least
squares deconvolution model. Let `Y_M` denote the normalized spatial expression
matrix restricted to marker genes, and let `S_M` denote the marker-gene by
cell-type reference matrix. For each spot, PASSAGE estimates a nonnegative
cell-type weight vector whose elements sum to one. In matrix form, the fitted
weights solve:

```text
min_W || W S_M - Y_M ||_F^2
subject to W_ik >= 0 and sum_k W_ik = 1 for every spot i.
```

The optimization was implemented using projected gradient descent with
row-wise projection onto the probability simplex. The default maximum number of
iterations was 350. The resulting spot-level cell-type proportion matrix was
used as a covariate matrix in the H2 and H3 analyses.

## Background spatial factors

Background spatial factors were estimated to capture broad tissue-level
transcriptional structure that is not specific to the pathway being tested. For
each Visium dataset, genes not belonging to the Hallmark collection were used as
candidate background genes. The most variable background genes were selected,
with a default maximum of 2,500 genes. Their log-normalized expression was
residualized against the estimated cell-type proportions, standardized, and
decomposed by singular value decomposition. The leading background factor scores
were retained as broad spatial expression covariates. The default number of
background factors was 6.

These background factors were used in H3 to make the pathway test more
stringent. They are intended to adjust for broad tissue architecture and global
spatial expression patterns, while leaving pathway-specific residual spatial
structure to be evaluated by PASSAGE.

## Covariate designs

PASSAGE was evaluated under nested covariate designs. The baseline technical
and spatial covariate matrix included:

- log total UMI count,
- log detected gene count,
- scaled spatial x-coordinate,
- scaled spatial y-coordinate,
- scaled squared x-coordinate,
- scaled squared y-coordinate,
- scaled x-by-y interaction.

An intercept is added internally by the design preparation functions where
required.

The nested designs were:

```text
H1 design: technical and low-order spatial covariates.
H2 design: H1 design plus Wu deconvolved cell-type proportions.
H3 design: H2 design plus background spatial factors.
H4 design: H2 design plus a two-level region indicator for region testing.
```

The H4 design is used with a separate region-enrichment test. The default region
definition was an epithelial-enrichment median split derived from the Wu
deconvolution output. Specifically, epithelial enrichment was defined from the
sum of malignant epithelial and normal epithelial proportions when available.
Alternative two-level region definitions can be substituted, including
pathology-derived or image-derived annotations.

## PASSAGE spatial factor model

For a pathway containing `p` genes measured across `n` spatial spots, PASSAGE
models the normalized expression matrix as:

```text
Y = X B + U A' + E.
```

Here, `Y` is the `n x p` expression matrix, `X` is the `n x q` covariate matrix,
`B` is the `q x p` coefficient matrix for observed covariates, `U` is an
`n x K` matrix of latent spatial factor scores, `A` is a `p x K` matrix of
gene loadings, and `E` is residual expression. The columns of `U` represent
latent spatial processes evaluated at the observed spot coordinates. The
loading matrix `A` determines how strongly each pathway gene contributes to
each latent spatial process.

The model is low-rank in the gene dimension. This is motivated by the
expectation that pathway-level spatial activity can often be represented by a
small number of coordinated spatial programs rather than one independent
spatial process per gene. Gene-specific residual variances and residual
covariance summaries are estimated after fitting the latent spatial component.

## Spatial covariance model

Each latent spatial factor is modeled using a spatial covariance kernel
evaluated at the observed spot coordinates. The implementation supports the
following kernels:

- Matern 1/2,
- Matern 3/2,
- Matern 5/2,
- exponential,
- Gaussian.

The production workflow used the Matern 3/2 kernel as the default. Spatial
range parameters were estimated over a grid determined from the tissue
coordinates. A Vecchia approximation was used for efficient spatial likelihood
and score calculations, with a default of 20 nearest neighbors.

## Rank selection

The number of latent factors can be fixed or selected dynamically. The dynamic
rank rule selects the smallest rank whose cumulative residual singular-value
variance exceeds a specified threshold, subject to a maximum rank. The
production all-engine workflow used:

```text
rank selection method = variance threshold
variance threshold = 0.90
maximum rank = 12
```

This allows the model rank to adapt to pathway complexity while limiting the
effective degrees of freedom of the latent factor representation.

## Spatial factor fitting engines

Several fitting engines were implemented to estimate the latent spatial factor
model. These engines return a common fitted object where possible, allowing the
same downstream tests and metrics to be applied across engines.

### PCA plus Vecchia GP engine

The PCA engine first residualizes pathway expression with respect to the chosen
covariate matrix. Singular value decomposition is then applied to the residual
expression matrix. The leading residual principal component score vectors are
treated as latent spatial factors, and spatial covariance parameters are fitted
for each factor using the Vecchia approximation. Gene loadings and residual
variances are refit after factor estimation.

This engine is computationally efficient and provides the default PASSAGE
factorization. It differs from simple PC adjustment because the PCs are used as
estimated pathway spatial factors and are subsequently assigned spatial
covariance structure.

### Spatial-basis ridge engine

The spatial-basis engine constructs radial-basis and polynomial features from
the spatial coordinates. The basis is residualized against the covariate design,
and a ridge-smoothed multigene expression surface is fitted. The smoothed
surface is then compressed by singular value decomposition to obtain latent
spatial factors and gene loadings.

This engine estimates smooth spatial expression structure directly from
coordinates before performing low-rank compression.

### Smoothed PCA engine

The smoothed PCA engine starts with residual PCA factors and then smooths the
factor scores over spatial coordinates using a ridge-penalized spatial basis.
Gene loadings are refit after smoothing.

This engine retains the dominant residual expression directions while reducing
non-spatial noise in the factor scores.

### NMF engine

The nonnegative matrix factorization engine shifts the residual expression
matrix to a nonnegative scale and fits a parts-based factorization. The fitted
NMF scores can be spatially smoothed, and signed gene loadings are refit on the
original residual expression scale.

This engine provides a non-orthogonal, parts-based alternative to PCA-like
factorizations.

### Alternating GP engine

The alternating GP engine alternates between updating latent factor scores and
gene loadings while applying a spatial covariance penalty to the factor scores.
This engine approximates direct spatial factor fitting while remaining more
computationally tractable than a full joint likelihood model.

### Joint TMB spatial factor engine

A joint likelihood engine was also implemented using Template Model Builder.
This engine directly fits a latent spatial factor model with Vecchia spatial
priors, gene-specific residual variances, and a lower-triangular loading
parameterization for identifiability. The TMB engine supports nested model
comparisons and model-implied pathway spatial variance summaries.

Because the downstream competitive and coherence testing infrastructure was
designed around the common fitted PASSAGE engine object, the production
all-engine workflow emphasized the PCA, spatial-basis, smoothed-PCA, NMF, and
alternating-GP engines. The TMB engine provides a direct likelihood-based
extension for future sensitivity analyses.

## Self-contained PASSAGE score test

For a pathway and a selected covariate design, PASSAGE tests the null hypothesis
that the pathway has no residual spatial component:

```text
H0: U A' = 0 after adjustment for X.
```

Equivalently, the test asks whether the genes in the pathway carry spatially
structured residual expression after conditioning on the selected covariates.
The score statistic is computed from residualized pathway expression and the
fitted factor-specific spatial kernels. The implementation evaluates multiple
factor-weighting schemes, including equal, variance-based, and range-based
weights, and reports an omnibus pathway p-value.

Permutation calibration was used for production pathway analyses. For each
pathway, the default production configuration used 999 permutations. A
moment-based p-value is also available as a faster approximation and was used in
some diagnostic and gene-level steps.

## Nested pathway hypotheses

PASSAGE implements nested self-contained pathway tests corresponding to H1, H2,
and H3.

The H1 null is:

```text
H1 H0: the pathway has no residual spatial signal after adjustment for
technical covariates and low-order coordinate effects.
```

The H2 null is:

```text
H2 H0: the pathway has no residual spatial signal after adjustment for
technical covariates, low-order coordinate effects, and estimated cell-type
composition.
```

The H3 null is:

```text
H3 H0: the pathway has no residual spatial signal after adjustment for
technical covariates, low-order coordinate effects, estimated cell-type
composition, and broad background spatial factors.
```

H1 is the most sensitive self-contained test. H2 asks whether pathway spatial
signal remains after accounting for cell composition. H3 is a stricter
conditional test that asks whether pathway spatial signal remains after
additional adjustment for broad background spatial structure.

The implementation also reports decomposition summaries comparing H1, H2, and
H3 statistics. These summaries quantify the share of spatial signal associated
with cell-type adjustment, background-factor adjustment, and pathway-specific
residual signal under the fitted model.

## Gene-level spatial testing and ACAT aggregation

As a gene-level baseline, PASSAGE computes a spatial p-value for each gene
appearing in at least one tested pathway. Each gene is tested using the same
spatial score-testing framework with H1 covariate adjustment. For a pathway,
gene-level spatial p-values are aggregated using the Cauchy combination test
(ACAT).

The SVG-ACAT baseline tests whether a pathway contains one or more genes with
strong gene-level spatial signal. This baseline is useful for comparison with
gene-centric spatial variability workflows, but it does not require pathway
genes to form a coordinated spatial program.

## Pathway spatial variance metrics

PASSAGE reports several pathway-level spatial variance and coherence metrics.
These metrics are intended to complement p-values and to support pathway
ranking by effect size.

### Mean proportion spatial variance

For each pathway gene, the fitted spatial variance is compared with the total
fitted variance. The mean proportion spatial variance is the average of this
quantity across genes:

```text
mean_propSV = mean_j spatial_variance_j /
              (spatial_variance_j + residual_variance_j).
```

This metric is simple and interpretable, but it is gene-marginal. It does not
fully account for covariance among pathway genes.

### Covariance-weighted spatial variance

Let `S_spatial` denote the covariance matrix of the fitted spatial pathway
signal and let `S_residual` denote a shrinkage estimate of the residual
expression covariance matrix. PASSAGE computes trace-based spatial variance:

```text
cwPVE_trace = trace(S_spatial) /
              trace(S_spatial + S_residual).
```

It also computes generalized spatial variance summaries from the covariance
structure of `S_spatial` relative to `S_spatial + S_residual`. These summaries
ask how much of the multivariate pathway covariance is attributable to fitted
spatial structure.

### Effective pathway spatial variability

Effective pathway spatial variability summarizes the contribution of latent
spatial factors while incorporating spatial scale. Factors with stronger
loadings and more spatially effective range parameters contribute more to the
metric than weak or very local factors.

### Coherence-weighted effective pathway spatial variability

Coherence-weighted effective pathway spatial variability, abbreviated cEPSV,
extends effective pathway spatial variability by weighting each spatial factor
by pathway loading coherence. A factor receives more weight when many genes in
the pathway load coherently on the same spatial factor and less weight when the
factor is driven by isolated genes.

Conceptually:

```text
cEPSV = sum_k spatial_contribution_k
              x spatial_scale_weight_k
              x pathway_loading_coherence_k
        / total_pathway_variation.
```

cEPSV was implemented to distinguish coordinated pathway-level spatial programs
from gene sets that contain individual spatially variable genes without strong
shared pathway structure.

## Competitive pathway testing

The Conditional Competitive Null evaluates whether a pathway's spatial metric
is larger than expected for matched background gene sets:

```text
H0: after conditioning on the selected covariates, the pathway is no more
spatially organized than expression-matched background gene sets of the same
size.
```

Background genes are matched to pathway genes using expression and detection
features. Matched random gene sets are sampled repeatedly, and the selected
pathway metric is recomputed for each sampled set. The empirical competitive
p-value is computed as:

```text
p = (1 + number of null statistics at least as large as observed) /
    (1 + number of permutations).
```

The competitive testing framework supports multiple metrics, including
mean proportion spatial variance, effective pathway spatial variability, cEPSV,
and covariance-weighted PVE summaries. Production competitive analyses used 999
matched-set permutations where configured.

## Pathway coherence testing

The Pathway Coherence Null evaluates whether the genes in a pathway are more
spatially coordinated than expected for matched background gene sets:

```text
H0: pathway genes are not more spatially coherent than matched random gene
sets.
```

Implemented coherence statistics include the spatial fraction aligned with the
first pathway expression axis and the mean absolute correlation of the fitted
spatial signal among pathway genes. As with competitive testing, significance is
assessed by matched-set permutation.

This test is designed to distinguish a coherent spatial pathway program from a
pathway that is significant only because it contains several unrelated
spatially variable genes.

## Region-enrichment testing

The Region-Enrichment Null evaluates whether fitted pathway spatial signal is
preferentially enriched across a specified two-level tissue-region contrast:

```text
H0: the pathway's fitted spatial signal is not differentially enriched across
the two region labels.
```

For each latent factor, PASSAGE computes region-specific spatial quadratic
forms using residualized pathway expression and the fitted spatial kernel. The
test statistic is a weighted sum of absolute differences between the two
region-specific quantities across latent factors. P-values are computed by
permuting region labels.

The default region contrast in the breast cancer workflow was an exploratory
epithelial-enrichment median split derived from deconvolved epithelial
proportions. The same method can be applied to histology-derived, manual, or
image-derived region labels.

## Multiple testing

For pathway-level analyses, p-values were adjusted within each dataset,
pathway collection, engine, and hypothesis using the Benjamini-Hochberg false
discovery rate procedure. Separate adjusted p-values were computed for H1, H2,
H3, H4, and SVG-ACAT outputs. Competitive and coherence p-values can be
adjusted analogously within each analysis stratum.

## Computational implementation

The method was implemented in R as part of the PASSAGE codebase. Core
functions include:

```text
passage_fit_engine_pca()
passage_fit_factor_engine()
passage_score_test()
passage_pathway_covariance_metrics()
passage_conditional_competitive_test()
passage_pathway_coherence_test()
passage_region_enrichment_test()
passage_fit_joint_factor_tmb()
```

The main 10x MSigDB workflow script is:

```text
scripts/run_10x_msigdb_conditional_passage.R
```

The deconvolution workflow is:

```text
scripts/deconvolve_10x_with_wu_reference.R
```

The Biowulf production workflow was configured to run Hallmark and KEGG
pathways across the two breast cancer Visium datasets and across the compatible
PASSAGE fitting engines:

- PCA plus Vecchia GP,
- spatial-basis ridge,
- smoothed PCA,
- NMF,
- alternating GP.

The production configuration used 999 permutations for pathway score tests and
999 matched-set permutations for competitive testing where configured.

## Simulation and software validation

Simulation scripts were implemented to evaluate null calibration and runtime
under controlled settings. These simulations include independent and correlated
null expression structures and permutation-based calibration with 999
permutations. The software test suite includes smoke tests for the PASSAGE
workflow, tests for the fitting engines, tests for competitive and coherence
metrics, and tests for the TMB joint factor implementation.

No simulation or real-data results are reported in this section. Simulation
outputs and empirical breast cancer pathway results should be summarized in a
separate Results section after the final analysis set is fixed.
