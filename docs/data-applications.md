# Breast and Kidney Cancer Applications

These workflows reproduce the packaged breast cancer and kidney cancer
applications on NIH Biowulf. They download public spatial data, construct the
six standard PASSAGE inputs, run self-contained screens and competitive H3,
and aggregate pathway and driver results. All computation is submitted to
Slurm; do not execute the R analysis programs on a login node.

## Requirements

- a Biowulf account with `sbatch`, R, Python, and MAGMA modules available;
- a PASSAGE checkout on Biowulf;
- space for public expression matrices, references, and results;
- internet access from the download jobs;
- `msigdbr` and the R packages used by the input-preparation scripts.

Create or update the code checkout:

```bash
CODE=/data/Dutta_lab/SPATH/PASSAGE
git clone https://github.com/diptavo/PASSAGE.git "$CODE"
# For an existing checkout:
git -C "$CODE" pull --ff-only
```

The launchers print every Slurm job ID. Monitor them with `squeue -u "$USER"`
and inspect each application's `logs/` directory if a dependency does not run.

## Breast Cancer

The breast application uses the 10x Genomics FFPE breast cancer Visium sample
and two broad cell-type covariate panels (`tumor_tme_broad` and
`normal_tissue_broad`). Copy the workflow into a results root and submit it:

```bash
CODE=/data/Dutta_lab/SPATH/PASSAGE
BREAST_ROOT=/data/Dutta_lab/SPATH/PASSAGE_breast_application
mkdir -p "$BREAST_ROOT"
rsync -a "$CODE/workflows/biowulf/biowulf_cancer_panel/" "$BREAST_ROOT/"
bash "$BREAST_ROOT/submit_pipeline.sh" \
  "$BREAST_ROOT" "$CODE" breast 999 9999
```

The final two numbers are the self-contained permutation count and competitive
Monte Carlo count. Use `199 999` for a smoke run. For the full four-cancer
panel, replace `breast` with `breast,cervical,prostate,lung`.

Primary breast outputs are:

```text
results/passage_cancer_summary/all_sample_competitive_h3.csv
results/passage_cancer_summary/competitive_h3_pathway_summary.csv
results/passage_cancer_summary/summary.md
```

## Kidney Cancer

The kidney application uses three RCC spatial samples (KC1, KC2, and KC3) and
two broad kidney cell-type covariate panels per sample. Submit the core spatial
analysis with:

```bash
CODE=/data/Dutta_lab/SPATH/PASSAGE
KIDNEY_ROOT=/data/Dutta_lab/SPATH/PASSAGE_kidney_application
mkdir -p "$KIDNEY_ROOT"
rsync -a "$CODE/workflows/biowulf/biowulf_kidney_rcc_gwas/" "$KIDNEY_ROOT/"
export PASSAGE_GWAS_GLOB='/data/Renal_GWAS_2022_exp/sumstats/meta_20230322/gwas_ssf/meta/meta.multianc.*.sumstats.tsv.bgz'
bash "$KIDNEY_ROOT/submit_pipeline.sh" \
  "$KIDNEY_ROOT" "$CODE" 999 9999
```

`PASSAGE_GWAS_GLOB` is recorded by the download job for the optional validation
stage. Primary kidney outputs are:

```text
results/passage_kidney_summary/all_sample_competitive_h3.csv
results/passage_kidney_summary/competitive_h3_pathway_summary.csv
results/passage_kidney_summary/summary.md
```

After the core pipeline completes successfully, submit post-hoc RCC, clear-cell
RCC, and papillary RCC GWAS validation:

```bash
bash "$KIDNEY_ROOT/submit_gwas_validation.sh" \
  "$KIDNEY_ROOT" \
  /data/Dutta_lab/REF/EUR \
  /data/Dutta_lab/tools/NCBI38.gene.loc
```

This stage formats the three summary-statistic files, performs MAGMA gene
analysis, tests PASSAGE pathway and driver gene sets, and writes results under
`results/magma_passage_sets_summary/`.

## Cell-Type Adjustment

The current application inputs use two marker-derived broad cell-type score
panels for each sample. The cancer-panel download job also stages raw scRNA-seq
reference archives, but the production scripts do not yet fit a probabilistic
cell-type deconvolution model from those archives. Thus, `Z_CT` controls broad
composition signatures; it should not be described as reference-based cell
fractions without replacing it with validated deconvolution estimates.

To use external proportions, retain the spot order in each input bundle,
replace `Z_CT` with the aligned proportion matrix, and rerun stages 02 through
04. Rows must match `Y` and `coords` exactly.

## Reading Results

Use `score_z` as the primary competitive statistic and
`competitive_score_fdr_global` as the study-wide multiplicity adjustment.
`score_robust_z` is a sensitivity analysis. cEPSV, coherence, effective pathway
size, and spatial-variance columns are descriptive biological summaries rather
than confirmatory p-value statistics.

```r
x <- read.csv("results/passage_kidney_summary/all_sample_competitive_h3.csv")
primary <- subset(x, statistic == "score_z")
primary <- primary[order(primary$competitive_score_fdr_global), ]
head(primary)
```

The empirical competitive p-value cannot be smaller than `1 / (B + 1)`.
Generalized-Pareto tail values should be used only when the corresponding tail
fit diagnostics pass and must be reported alongside the empirical p-value.
