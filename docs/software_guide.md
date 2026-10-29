# Software guide

PASSAGE provides R functions and two command-line programs. Start with the
pathway workflow for self-contained testing, or use the competitive workflow
when the question is enrichment relative to matched genes. The
[methods guide](method_choices.md) explains the distinction.

## Install

```r
install.packages("remotes")  # if needed
remotes::install_github("diptavo/PASSAGE")
library(PASSAGE)
packageVersion("PASSAGE")
```

For a source archive, use
`install.packages("PASSAGE_0.5.0.tar.gz", repos = NULL, type = "source")`.
The core package requires R >= 4.1 and `Matrix`. Install `irlba` for fast global
PCA with large backgrounds. Descriptive GP fitting needs `GpGp` and `fields`;
the optional compiled joint-factor engine needs `TMB`, `RcppEigen` and a working
C++ compiler. Record the package version and `sessionInfo()` with each analysis.

## Kidney example: estimation and both testing questions

The package includes a prepared public KC2 kidney cancer sample with **2,166
spots and 17,630 genes**. Every prepared gene is retained so that competitive
tests have a background beyond the target pathways. Coordinates, normalized
expression, raw counts, QC measurements, reference-informed composition
estimates and the adjustment design are already aligned.

Three Hallmark sets are included: hypoxia (190 observed genes), angiogenesis
(32), and inflammatory response (197). Their complete gene lists are retained,
including genes absent from this expression matrix. See the
[data provenance and licenses](../inst/extdata/kidney_KC2_SOURCE.md).

### Run from a terminal

Download the repository once and install the package:

```bash
git clone https://github.com/diptavo/PASSAGE.git
cd PASSAGE
Rscript -e 'install.packages(c("remotes", "irlba"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_local(".", upgrade="never")'
```

If the repository is already on your computer, enter its `PASSAGE` directory
and update it with `git pull --ff-only` instead of cloning again. Run all three
pathways through estimation and testing:

```bash
Rscript inst/examples/run_kidney.R
```

This uses three pathway PCs and 199 null draws and writes `results/kidney_all/`.
The command performs:

1. Covariate-adjusted PCA within each pathway, saving scores and loadings.
2. A self-contained GP covariance test using those scores, with a shared
   64-feature kernel cache across pathways.
3. Fixed-kernel trace and largest-root tests, combined into one corrected
   omnibus p-value per pathway, as a sensitivity analysis.
4. Global PCA across all 17,630 genes and competitive testing against unique
   matched gene sets, excluding the target pathway from its own background.

For more Monte Carlo resolution, use a new output directory:

```bash
Rscript inst/examples/run_kidney.R all 9999 results/kidney_9999
```

The `all` command includes every route above. To run only one testing route:

```bash
Rscript inst/examples/run_kidney.R competitive 9999 results/kidney_competitive
Rscript inst/examples/run_kidney.R gp 9999 results/kidney_gp
Rscript inst/examples/run_kidney.R fixed 9999 results/kidney_fixed
```

The optional final arguments choose the pathway estimation method and rank:

```bash
Rscript inst/examples/run_kidney.R gp 199 results/kidney_fa factor_analysis 3
Rscript inst/examples/run_kidney.R gp 199 results/kidney_pca95 pca auto
```

`pca auto` retains 95% of adjusted pathway variance. The example's default
`pca 3` is a smaller, prespecified model for a quick run; it tests only the
retained three-component representation. It does not change the core function's
95% default. Factor analysis needs a fixed rank. These choices affect the
pathway-specific route; competitive analysis uses its own global six-factor PCA
with 20 Vecchia neighbors. The example specifies ranges of 0.05, 0.1, 0.2 and
0.4 for the GP and competitive routes, and 0.25 for the fixed-kernel route.

### Run from R without a repository checkout

After installing PASSAGE and `irlba`, the same example is available through the
installed package:

```r
library(PASSAGE)
source(system.file("examples", "run_kidney.R", package = "PASSAGE"))
run_kidney()
```

Choose a new output directory for another run. The R entry point also exposes
the kernel feature count and optional descriptive GpGp fits:

```r
run_kidney(mode = "gp", B = 199, n_features = 128,
           output = "results/kidney_gp128")

# Optional descriptive fits; the GP test itself does not need these packages.
install.packages(c("GpGp", "fields"))
run_kidney(mode = "gp", B = 199, descriptive_gp = TRUE,
           output = "results/kidney_descriptive")
```

To inspect the data or use the individual functions below:

```r
b <- passage_validate_bundle(system.file("extdata", "kidney_KC2.rds", package = "PASSAGE"))
dim(b$Y)
names(b$pathway_sets)
colnames(b$X)
```

### Read the results

`summary.csv` has one row per pathway and test: nine rows for `all`. It records
the raw p-value, BH/BY adjustments across the three pathways within each test,
and BH/BY adjustments across all nine reported pathway-test hypotheses. Select
the testing family before interpreting results. The GP test is the example's
primary self-contained test; the fixed-kernel result is a sensitivity analysis.
Competitive testing answers a different question.

Each pathway directory contains the fitted components, score/loadings tables,
and test objects. `competitive_summary.csv` contains competitive p-values,
effect-size and coherence summaries, and candidate driver genes;
`competitive_result.rds` also retains matched-null results. The top-level
configuration records the input checksum, design, coordinate scaling, settings
and runtime. `sessionInfo.txt` records the R environment, and `analysis.complete`
is written only after every requested step succeeds. Existing nonempty output
directories are rejected to protect earlier results.

The 199-draw run checks execution; it has a raw Monte Carlo p-value floor of
0.005. Increasing to 9,999 draws improves that floor to 0.0001. The fixed-test
omnibus has twice that floor because it corrects over two statistics. More
draws do not fix model assumptions.
The Gaussian-reference tests remain sensitive to spatially varying noise
variance. Composition covariates here were estimated from the same counts and
are treated as fixed, without propagating their uncertainty. Review QC,
calibration and replication before making biological claims. This example
covers spatial estimation and testing; GWAS integration still requires separate
GWAS statistics, gene mapping and appropriate LD information.

## Prepare an input bundle

The bundle is an ordinary R list saved as an RDS file:

```r
bundle <- list(
  Y = expression_matrix,       # spots x genes, normalized expression
  coords = coordinate_matrix, # spots x spatial dimensions
  X = design_matrix,          # spots x adjustment covariates; NULL is intercept-only
  pathway_sets = list(
    pathway_A = paste0("GENE", 1:20),
    pathway_B = paste0("GENE", 21:40)
  ),
  sample = "sample_01"
)
bundle <- passage_validate_bundle(bundle)
saveRDS(bundle, "sample.rds")
```

The named objects above stand for your prepared data. Use unique gene names in
`colnames(Y)` and unique spot IDs in `rownames(Y)`. The row names of `coords`
and a supplied `X` must match those IDs in exactly the same order. Values must
be finite, coordinates unique and the design full rank. Validation checks
these structural requirements; it does not normalize counts, map gene IDs or
establish the statistical assumptions.

An intercept is added if absent. If composition fractions sum to one, omit one
fraction from the design when using an intercept to avoid collinearity. An
optional `cell_type_proportions` matrix records all aligned fractions and is
checked for valid row sums, but does not automatically add adjustment columns
to `X`.

The validator requires at least 10 spots and five genes in the bundle. Pathway
fitting requires at least two observed genes with nonzero adjusted variance.
Competitive analysis defaults to observed pathway sizes of 15–500 genes and
needs a sufficiently large matching background.

## Estimate and test one pathway

```r
library(PASSAGE)
b <- passage_validate_bundle("sample.rds")
pathway <- "pathway_A"

# Use a common scale across axes to retain the tissue's aspect ratio.
xy <- sweep(b$coords, 2, apply(b$coords, 2, min), "-")
xy <- xy / max(xy)

fit <- passage_fit_pathway(
  b$Y, b$pathway_sets[[pathway]], X = b$X,
  method = "pca", variance_threshold = 0.95
)

test <- passage_test_pathway_gp(
  fit, coords = xy,
  ranges = c(0.05, 0.1, 0.2, 0.4),
  kernel_method = "lowrank", n_features = 64,
  B = 9999, seed = 1
)

fit$K
fit$missing_genes
fit$dropped_genes
test$p
test$component_fits
saveRDS(list(fit = fit, test = test), "pathway_result.rds")
```

The R API retains the coordinate units supplied to it. The example scales the
longest axis to span one; the range grid is in those units. The command-line
program performs and records this scaling automatically.

Use `K = 2` to fix the PCA rank. For factor analysis, set
`method = "factor_analysis", K = 2`. FA needs an identifiable factor count and
a nonsingular gene correlation matrix; a rank chosen for PCA is not
automatically suitable for FA.

`kernel_method = "dense"` avoids the kernel approximation and is preferable
for small datasets where it is affordable. Its default limit is 1,000 spots.
The low-rank option retains all spots and uses `n_features` landmarks per range.
Assess approximation sensitivity before committing to a large analysis.

The default spatial-fraction grid is
`c(0, .01, .025, .05, .1, .25, .5, .75, .9, .97)`. A custom grid must include
zero and contain values below one. Range and fraction grids are searched for
the observed scores and every null draw. Choose them before inspecting results.

### Fixed-kernel alternative

To test alignment with a prespecified kernel rather than fit its range, use:

```r
fixed <- passage_test_pathway_fit(
  fit, coords = xy, kernel_range = 0.25, B = 9999, seed = 1
)
fixed$p                       # trace and largest-root p-values
fixed$p_omnibus_bonferroni     # correction over those two tests only
```

This dense reference has a default limit of 3,000 spots. Its omnibus correction
does not adjust for testing multiple pathways or samples. Designate a primary
test in advance if running both fixed-kernel and GP analyses.

### Reuse kernels across pathways

```r
cache <- passage_prepare_gp_kernel(
  xy, X = b$X, ranges = c(0.05, 0.1, 0.2, 0.4), n_features = 64
)
test <- passage_test_pathway_gp(fit, kernel_cache = cache, B = 9999, seed = 1)
```

Reuse the cache only with fits from the same spots, in the same order, with the
same design. Cache settings govern the test; conflicting explicit settings
cause an error.

### Optional descriptive GP fits

```r
install.packages(c("GpGp", "fields"))  # once, if this option is needed
descriptive <- passage_fit_pathway_gp(fit, xy)
descriptive$summary
```

This reports GpGp covariance estimates and convergence. It neither supplies
inputs to the tests above nor returns a spatial-test p-value. The GP comparison
already fits its own covariance parameters and does not require this step.

## Competitive pathway analysis

```r
b <- passage_validate_bundle("sample.rds")
competitive <- passage_analyze_bundle(b, B = 9999, seed = 1)
competitive$summary
competitive$coverage
saveRDS(competitive, "competitive_result.rds")
```

The bundle must include the background gene universe, not just the tested
pathway. This route uses a global PCA engine and unique matched comparison
sets. The main result fields are `p_competitive`, `q_BH`, `q_BY`, `score_excess`
and `candidate_drivers`. Adjustments cover eligible pathways in that call; if
the family includes several samples, collect and adjust the relevant p-values
together.

If matching fails, inspect background coverage and the strata before changing
the matching design. To supply bins, pass a `gene_bins` data frame with `gene`
and `bin`, plus an optional `module`. Its rows must exactly follow the
nonconstant expression columns. Do not silently include the tested genes in
their own background.

## Command line

From the repository's `PASSAGE` directory after installing the package, run a
pathway GP comparison using a pathway name present in your bundle:

```bash
Rscript inst/scripts/passage_pathway_cli.R sample.rds pathway_A gp_results \
  9999 1 pca auto 0.25 false gp lowrank 64
```

For the fixed-kernel route:

```bash
Rscript inst/scripts/passage_pathway_cli.R sample.rds pathway_A fixed_results \
  9999 1 pca auto 0.25 false subspace
```

Each pathway output directory must be new or empty. The argument order is:

```text
BUNDLE PATHWAY OUTPUT B SEED METHOD K RANGE DESCRIPTIVE_GP
TEST GP_KERNEL FEATURES RANGES FRACTIONS
```

| Argument | Meaning |
|---|---|
| `B`, `SEED` | Number of null draws and random seed |
| `METHOD`, `K` | `pca auto` for 95% variance retention, `pca 2` for fixed rank, or `factor_analysis 2` |
| `RANGE` | Fixed-kernel range; retained as a positional argument but unused by the GP test |
| `DESCRIPTIVE_GP` | `true` additionally fits the separate GpGp models; otherwise `false` |
| `TEST` | `subspace` or `gp`; the default is `subspace` |
| `GP_KERNEL`, `FEATURES` | `lowrank` or `dense`, and the low-rank feature count |
| `RANGES`, `FRACTIONS` | Optional comma-separated grids for the GP test |

The final two arguments default to the grids described above. The older
competitive CLI takes a bundle and a separate RDS containing a named list of
gene sets:

```bash
Rscript inst/scripts/passage_cli.R sample.rds pathways.rds competitive_results 9999 1
```

Use a fresh directory for competitive results as well; that older CLI does not
reject an existing directory.

## Outputs and interpretation

| Pathway CLI output | Contents |
|---|---|
| `pathway_fit.rds` | Estimation object, selected rank and gene coverage |
| `component_scores.csv`, `component_loadings.csv` | Spot activity and gene contributions |
| `pathway_test.rds` | Statistic, null draws, p-values and test settings |
| `pathway_summary.csv` | `p_gp` for the GP route, or `p_trace`, `p_root`, `p_omnibus` for the fixed route |
| `gp_component_fits.csv` | GP comparison estimates and grid-boundary flags, when applicable |
| `configuration.rds`, `sessionInfo.txt` | Input checksum, coordinate scaling, software version and settings |
| `analysis.complete` | Written after successful execution |

The competitive CLI writes `passage_result.rds`, `pathways.csv`,
`pathway_coverage.csv`, `sessionInfo.txt` and `analysis.complete`.

For a quick execution check, use `B = 39` or `199`. For inference, the minimum
attainable p-value is `1/(B+1)`: 9,999 draws allow 0.0001. Choose enough draws
for the multiple-testing threshold rather than treating that floor as a precise
tail probability. The single-pathway functions do not adjust across pathways.

Successful execution does not establish calibration. In particular, the current
Gaussian reference is not robust to spatially varying noise variance; see the
[methods guide](method_choices.md#validation-and-current-limits). Candidate
drivers, component loadings and fitted spatial parameters need biological
interpretation and replication. GWAS integration is a separate analysis.

## Package checks

From a source checkout:

```bash
R CMD build .
R CMD check --no-manual PASSAGE_0.5.0.tar.gz
```

The test suite covers factor reconstruction, matched-set sampling, input
alignment, pathway estimation, GP likelihood calculations and null rotations.
Optional engines are exercised when their dependencies are available.
