# Preparation, attribution and data licenses

These bundles are derived from public spatial-transcriptomics datasets. They
retain the complete prepared expression matrices frozen on October 22, 2026.
The public export on November 5, 2026 removes duplicate arrays and local
filesystem paths, while retaining the analysis arrays and available metadata.
There is no new spot filtering, gene filtering, normalization or composition
estimation in this export.

## Public expression and annotation sources

**Kidney KC1–KC3 and lung LC1–LC5.** Dawo, Sebastian; Nonchev, Kalin; Silina,
Karina (2025). *10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung
(5) Cancer with Tertiary Lymphoid Structures*, version 1, January 9, 2025.
[Zenodo record 14620362](https://zenodo.org/records/14620362),
DOI [10.5281/zenodo.14620362](https://doi.org/10.5281/zenodo.14620362).
Distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Counts and coordinates are from the Space Ranger outputs in that deposit.
Existing spot annotations were extracted from its `h5ad_preprocessed` files
and aligned to the retained count barcodes; labels and missing values are
preserved, with annotation-availability flags.

**Breast FFPE.** 10x Genomics,
[Human Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE)](https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0).
The prepared files use Space Ranger 1.3.0 sample
`Visium_FFPE_Human_Breast_Cancer`, with 2,518 spots. CC BY 4.0.
This specific source replaces an imprecise historical whole-transcriptome
landing-page link; the underlying prepared matrix is unchanged.

**Breast Block A.** 10x Genomics,
[Human Breast Cancer (Block A Section 1)](https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0),
Space Ranger 1.1.0 sample `V1_Breast_Cancer_Block_A_Section_1`. CC BY 4.0.

**Prostate.** 10x Genomics,
[Human Prostate Cancer, Adjacent Normal Section with IF Staining (FFPE)](https://www.10xgenomics.com/datasets/human-prostate-cancer-adjacent-normal-section-with-if-staining-ffpe-1-standard-1-3-0).
The prepared matrix uses the Space Ranger 2.0.0 download for
`Visium_FFPE_Human_Prostate_IF`; the source landing page describes the original
1.3.0 processing. CC BY 4.0. This section is adjacent-normal tissue from a
prostate cancer case, not an independent primary prostate tumor.

**Cervical.** 10x Genomics,
[Human Cervical Cancer (FFPE)](https://www.10xgenomics.com/datasets/human-cervical-cancer-1-standard),
Space Ranger 1.3.0 sample `Visium_FFPE_Human_Cervical_Cancer`. CC BY 4.0.

**NSCLC P4.** Behera V; Raman AS; Garassino MC; Kline JP.
*Conserved principles of spatial biology define tumor heterogeneity and
immunotherapy sensitivity*, [GSE292299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292299),
sample [GSM8855706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8855706).
This sample is an NSCLC brain biopsy. The export derives from the public
processed gene-count matrix and spatial coordinates. The GEO record does not
state a dataset-specific Creative Commons license, so none is assigned here;
retain the source attribution and consult the
[GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
The study did not deposit raw human sequencing reads because of privacy
concerns. This folder contains no sequencing reads.

Each expression dataset retains its source terms. PASSAGE's MIT license covers
the software, not a blanket relicensing of these third-party datasets. These
exports are adaptations and do not imply endorsement by the source creators.

## Expression preparation

All bundles align in-tissue spot barcodes with the source filtered count
matrices and use normalized expression:

```text
Y[spot, gene] = log1p(10000 * count[gene, spot] / original_spot_library_size)
```

Library sizes are calculated before gene filtering. Duplicate gene symbols
were collapsed during preparation. Genes with effectively zero normalized
variance were removed from the final bundles.

The four original production bundles, KC1–KC3 and breast FFPE, retained genes
detected in at least three spots. Their preparation did not apply an additional
UMI or detected-gene spot threshold beyond the source in-tissue selection.

The remaining nine bundles used these prespecified filters:

- at least 500 counts and 200 detected genes per spot;
- at most 30% mitochondrial counts and finite coordinates;
- genes detected in at least `max(3, ceiling(0.005 * retained_spots))` spots.

These differences in original preparation are preserved and should be considered
when comparing samples. FFPE probe panels may omit mitochondrial genes, making
the mitochondrial fraction uninformative for some samples. Raw counts in the
export have the same retained genes and spots as `Y`; they are not the complete
unfiltered source feature-barcode matrix. Histology images remain at the source.

## Adjustment covariates

Each `X` includes an intercept and the retained, standardized technical
covariates. Candidate technical variables include log library size, detected
genes, and mitochondrial fraction; the original production preparation also
considered ribosomal fraction. Constant or linearly redundant columns were
removed. Inspect `colnames(b$X)` for the exact design of each sample.

Kidney and breast additionally include reference-informed composition estimates.
The single-cell reference accessions are
[GSE224630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630) for kidney
and [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) for
breast. Annotated-cell signatures were restricted to curated markers, and
fractions were estimated by simplex-constrained quadratic programming with an
unmodeled-background component. One composition column was omitted from the
design to avoid collinearity. Full fractions and the omitted category are
retained when present in the prepared input. Reference expression matrices are
not redistributed.

The fractions were estimated from the same spatial counts and are treated as
fixed by the analysis. Their uncertainty is not propagated automatically.
Lung, cervical and prostate inputs have **technical-only adjustment** and must
not be described as composition-adjusted. Keeping `X` preserves the original
analysis question; changing it changes what variation is adjusted away.

## Pathway definitions

All 50 lists are from the
[MSigDB human Hallmark collection](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp),
frozen with the prepared inputs on October 22, 2026. The exact MSigDB release
number was not recorded. These included lists, rather than a current download,
define the reproducible testing family.

Copyright 2004–2025 Broad Institute, Inc., Massachusetts Institute of Technology,
and Regents of the University of California. Hallmark gene sets are distributed
under CC BY 4.0 according to the
[MSigDB license terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
Complete gene-set definitions are retained even when some listed genes are
absent from an individual expression matrix.
