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

No expression matrices, reference panels, GWAS summary statistics, MSigDB
files, or derived results are stored in this repository. Users are responsible
for obtaining data under the source providers' licenses and access rules.

Never run computational analyses on the Biowulf login node. Submit the `.sbatch`
entry points through SLURM and use job arrays for pathway, bootstrap, and null
replicates.
