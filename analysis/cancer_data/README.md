# Prepared cancer datasets

This folder contains **13 complete prepared spatial-transcriptomics samples**,
covering **44,736 spots** across five tissue groups. Each `.rds` file can be
passed directly to `PASSAGE::passage_validate_bundle()` and
`PASSAGE::passage_analyze_bundle()` with PASSAGE 0.6.1 or later.

All retained genes are included for competitive testing. These are the full
prepared samples, with no additional spot or gene subsampling for publication.
Each sample includes all **50 human Hallmark gene sets**. The frozen inputs were
prepared on October 22, 2026; this public export preserves their expression,
coordinates, adjustment designs, counts, and existing annotations.

## Samples

| Sample ID | Spots | Genes | Adjustment |
|---|---:|---:|---|
| `kidney__KC1` | 4,173 | 17,208 | Technical + composition |
| `kidney__KC2` | 2,166 | 17,630 | Technical + composition |
| `kidney__KC3` | 3,217 | 16,571 | Technical + composition |
| `lung__LC1` | 2,013 | 16,731 | Technical |
| `lung__LC2` | 4,372 | 17,385 | Technical |
| `lung__LC3` | 3,514 | 15,200 | Technical |
| `lung__LC4` | 1,551 | 13,652 | Technical |
| `lung__LC5` | 2,009 | 14,134 | Technical |
| `lung__GSM8855706_NSCLC_P4` | 9,175 | 17,986 | Technical |
| `breast__BlockA` | 3,798 | 18,038 | Technical + composition |
| `breast__Visium_FFPE_Human_Breast_Cancer` | 2,518 | 16,782 | Technical + composition |
| `prostate__Visium_FFPE_Human_Prostate_IF` | 3,460 | 14,469 | Technical |
| `cervical__Visium_FFPE_Human_Cervical_Cancer` | 2,770 | 16,636 | Technical |

The prostate sample is an **adjacent-normal section from a prostate cancer
case**, as described by its source. The NSCLC P4 sample is a **brain biopsy**.
The tissue-group names identify the study groups and should not be taken to
mean that every section is a primary tumor from that organ. The samples come
from different studies and assays, rather than a harmonized patient cohort.

## Download one or all samples

In R, download a sample and verify its checksum:

```r
source("https://raw.githubusercontent.com/diptavo/PASSAGE/main/analysis/cancer_data/download_cancer_data.R")
files <- download_cancer_data("kidney__KC2")
```

For the full collection (about 375 MiB compressed), omit the sample argument:

```r
files <- download_cancer_data()
```

Files are saved under `cancer_data/`, with a subfolder for each tissue group.
Existing files with matching checksums are reused. To select several samples:

```r
files <- download_cancer_data(c("kidney__KC1", "kidney__KC2", "kidney__KC3"))
```

Compressed file size understates analysis memory use. The largest NSCLC sample
used about 20 GiB during the full-dimension smoke test; plan for at least 32 GB
of RAM for that sample. KC2 is a smaller starting example.

If you cloned the repository, the files already exist under
`analysis/cancer_data/`; no download helper is needed. These files are excluded
from source archives and the installed R package. Use the download helper or
`git clone` to obtain them. The smaller package example in `inst/extdata/`
contains the same full KC2 expression matrix but only three pathway definitions.

## Run competitive tests

Install PASSAGE and its fast PCA dependency once:

```r
install.packages(c("remotes", "irlba"))
remotes::install_github("diptavo/PASSAGE@v0.6.1", upgrade = "never")
```

After downloading KC2 as above:

```r
library(PASSAGE)
b <- passage_validate_bundle("cancer_data/kidney/kidney__KC2.rds")
dim(b$Y)
colnames(b$X)

selected <- b$pathway_sets[c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)]

result <- passage_analyze_bundle(
  b, pathways = selected, B = 199, seed = 1,
  reference = "universe",
  range_grid = c(0.05, 0.1, 0.2, 0.4),
  statistics = c("score_z", "score_activity_hotspot",
                 "score_robust_z", "score_sparse_topk_z"),
  combine = "minp"
)
result$summary[, c("pathway", "n_genes", "p_competitive", "q_BH")]
result$tests
saveRDS(result, "kidney_competitive.rds")
```

For a repository checkout, change only the input path to
`analysis/cancer_data/kidney/kidney__KC2.rds`. Select a different sample by
changing that path. Keep the complete expression matrix when selecting
pathways; it supplies the competitive background. Omit `pathways = selected`
to analyze all 50 included Hallmark sets.

`B = 199` is a quick execution check. Use more draws, for example `B = 9999`,
for a finer p-value resolution. The primary summary reports the prespecified
minimum-p combination and BH adjustment across the selected pathways. The
full-universe reference permits overlap with the target and assumes uniform
gene-set membership within matching strata under its null. Matching expression
and detection characteristics does not by itself establish that assumption
for a curated, co-regulated pathway.

## Run pathway estimation and a self-contained test

Using the same `b` object:

```r
xy <- sweep(b$coords, 2, apply(b$coords, 2, min), "-")
xy <- xy / max(xy)

fit <- passage_fit_pathway(
  b$Y, b$pathway_sets$HALLMARK_HYPOXIA, X = b$X, K = 3, seed = 1
)
kernel <- passage_prepare_gp_kernel(xy, X = b$X, n_features = 64)
test <- passage_test_pathway_gp(fit, kernel_cache = kernel, B = 199, seed = 1)
test$p
```

This uses the Gaussian-reference GP test. See the
[method guide](https://diptavo.github.io/PASSAGE/methods/) and
[unequal-noise guide](../../docs/robust_calibration.md) for its assumptions and
the wild-bootstrap option. For pathway-local competitive fitting, use
`estimation = "pathway"` in the bundle interface and select `gp_aggregations`
instead of the global `statistics` panel.

## Contents of each bundle

| Field | Contents |
|---|---|
| `Y` | Normalized expression, spots by genes, ready for PASSAGE |
| `coords` | Aligned two-column image-pixel coordinates |
| `X` | Aligned, full-rank adjustment design, including an intercept |
| `raw_counts` | Sparse genes-by-spots counts for every retained gene and spot |
| `pathway_sets` | Named list of 50 complete Hallmark definitions |
| `spot_data`, `gene_data` | Retained spot and gene metadata |
| `cell_type_proportions`, `cell_type_reference` | Saved composition estimates and their reference, where available |
| `pathology_annotations` | Original annotation values aligned by barcode, for KC1–KC3 and LC1–LC5 |
| `sample_metadata`, `source_provenance` | Tissue context, source attribution, licenses and frozen-input checksum |

Existing QC thresholds, composition diagnostics and sensitivity covariates are
also retained where available. A missing composition or pathology field means
that these data were not prepared for that sample; no values are imputed.
Source pathology annotations can contain missing values, and their coordinate
columns may use a different image resolution. Use `coords` for PASSAGE and
barcodes to align annotation labels.

The raw count matrix contains the retained genes, so its column sums need not
equal the original normalization denominators. Use the provided normalized
`Y` for direct reproduction. Do not normalize it a second time.

## Provenance and verification

- [Manifest](manifest.csv): filenames, dimensions, file sizes, MD5 and SHA256 checksums.
- [Sources](sources.csv): original prepared-input checksums, public accessions and tissue context.
- [Preparation and attribution](SOURCES.md): normalization, QC, adjustment and source licenses.
- [Validation results](validation.csv): full-dimension loading, matching and execution checks with PASSAGE 0.6.1.

`export_prepared.R` repackages a frozen prepared input and verifies exact equality
of the main analysis arrays before and after export. It does not recreate the
original reference construction from raw sequencing data. `validate_prepared.R`
checks all 50 pathway definitions and runs a small competitive and self-contained
analysis without reducing the spot or gene dimensions. Those smoke tests check
software execution, not biological findings or statistical calibration.
