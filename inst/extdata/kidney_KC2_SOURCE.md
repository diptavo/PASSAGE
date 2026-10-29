# KC2 example data

`kidney_KC2.rds` contains the complete prepared KC2 expression matrix: 2,166
spots and 17,630 genes. The full background is retained for competitive testing.
The file also contains raw counts for these genes, coordinates, a 13-column
adjustment design, composition estimates, technical QC measurements, and three
Hallmark gene sets. No additional spot or gene subset was selected for this
example.

## Source and attribution

Dawo, Sebastian; Nonchev, Kalin; and Silina, Karina (2025).
*10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung (5) Cancer
with Tertiary Lymphoid Structures*, version 1.
[Zenodo, doi:10.5281/zenodo.14620362](https://zenodo.org/records/14620362).
The source is distributed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
This example uses the public KC2 Space Ranger expression and tissue coordinates.

The three gene lists are from the
[MSigDB human Hallmark collection](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp).
Copyright 2004-2025 Broad Institute, Inc., Massachusetts Institute of Technology,
and Regents of the University of California; distributed under CC BY 4.0 under
the [MSigDB license terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
The lists were frozen in the prepared input on October 22, 2026. Its exact
MSigDB release was not recorded; the included lists define this reproducible
example and should not be described as the latest release.

Reference-informed composition estimates were calculated using an annotated
kidney/renal cancer single-cell reference from
[GSE224630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630), restricted
to 81 overlapping curated marker genes. Reference expression data are not
redistributed in this file. The estimation method and marker names are recorded
in `cell_type_reference`.

These data retain their source licenses; the package's MIT license applies to
its code. The PASSAGE preprocessing and public export are adaptations and do
not imply endorsement by the original data creators.

## Preparation

The prepared input matches in-tissue barcodes to the filtered Space Ranger
matrix, collapses duplicate gene symbols, and retains genes detected in at
least three spots with nonzero normalized variance. Expression is
`log1p(10000 * counts / original_library_size)`, using library sizes calculated
before gene filtering. The public export preserves every value of `Y`, `coords`
and `X` in the frozen prepared input.

`X` contains an intercept, standardized log library size and detected-gene count,
and ten standardized composition columns. The pericyte fraction is omitted as
the reference category; constant technical columns are omitted. The full eleven
fractions, including an unmodeled-background component, are retained in
`cell_type_proportions`. Fractions are estimated from the same spatial counts
using the external reference. Treating them as fixed covariates does not account
for estimation uncertainty.

There was no extra UMI or detected-gene spot filter in this sample's preparation.
Review `spot_qc` before using it for a scientific study. Spatial coordinates are
image pixels; the example runner translates them and divides both axes by the
same maximum span. Histology images and manual tissue annotations are available
from the source record but are not included here.

| Pathway | Genes in the supplied set | Genes in the expression matrix |
|---|---:|---:|
| HALLMARK_HYPOXIA | 200 | 190 |
| HALLMARK_ANGIOGENESIS | 36 | 32 |
| HALLMARK_INFLAMMATORY_RESPONSE | 200 | 197 |

The export removes duplicate expression arrays, unrelated gene sets, absolute
filesystem paths and job metadata. It retains raw counts as a genes-by-spots
sparse matrix; `Y` is the normalized spots-by-genes matrix used by PASSAGE.

The export script is `inst/examples/prepare_kidney.R`. It accepts the frozen
prepared input and a new output path, checks the input checksum, and verifies
the saved object. It repackages that input; it does not reconstruct the external
reference or the original preprocessing from raw downloads.

- Frozen prepared-input MD5: `ba665dc5727ba41edef9c0dd97658955`.
- Public example MD5: `755944e9845d917b6bdccadecef3492f`.

The [software guide](https://github.com/diptavo/PASSAGE/blob/main/docs/software_guide.md#kidney-example-estimation-and-both-testing-questions)
gives the commands for analysis.
