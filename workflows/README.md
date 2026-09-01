# Reproducibility Workflows

`biowulf/` contains the SLURM pipelines used for PASSAGE development and
validation. They cover breast cancer, kidney cancer with post-hoc GWAS
validation, DLPFC, statistic calibration, pathway condensation, input-bundle
construction, and driver stability.

These scripts are analysis records and templates, not installed package code.
Their default `/data/...` paths identify the original NIH environment. Most
entry points accept the project root and reference paths as command-line
arguments; inspect each `.sbatch` file before submission and replace defaults
for another cluster.

The breast and kidney workflows include dependency-aware launchers:

```text
biowulf/biowulf_cancer_panel/submit_pipeline.sh
biowulf/biowulf_kidney_rcc_gwas/submit_pipeline.sh
biowulf/biowulf_kidney_rcc_gwas/submit_gwas_validation.sh
```

See [`../docs/data-applications.md`](../docs/data-applications.md) for exact
Biowulf setup, submission, output, and interpretation instructions.

No expression matrices, reference panels, GWAS summary statistics, MSigDB
files, or derived results are stored in this repository. Users are responsible
for obtaining data under the source providers' licenses and access rules.

Never run computational analyses on the Biowulf login node. Submit the `.sbatch`
entry points through SLURM and use job arrays for pathway, bootstrap, and null
replicates.
