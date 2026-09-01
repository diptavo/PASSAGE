# PASSAGE spatialDLPFC Biowulf Run

This run area downloads and prepares the LIBD spatialDLPFC dataset for PASSAGE
competitive H3 pathway testing.

## Dataset

- Spatial data: `spatialLIBD::fetch_data("spatialDLPFC_Visium")`
- Cell-type reference: `spatialLIBD::fetch_data("spatialDLPFC_snRNAseq")`
- Rationale: the same project provides Visium sections and DLPFC snRNA-seq,
  which is the correct background reference for estimating spot cell-type
  proportions.

## Analysis Goal

For each DLPFC section, run PASSAGE pathway tests under:

- H1: spatially varying pathway without cell-type adjustment.
- H3: residual spatial pathway variation after adjusting for technical
  covariates and estimated cell-type proportions.

The primary competitive H3 p-value should use `score_z` or `score_robust_z`
with empirical calibration. Secondary metrics can include pathway coherence,
activity hotspot strength, sparse/top-k driver concentration, and cEPSV as a
descriptive quantity only.

## Files

- `scripts/00_download_spatialDLPFC.R`: downloads Visium and snRNA-seq data.
- `scripts/00_download_spatialDLPFC.sbatch`: Biowulf batch wrapper.
- `scripts/01_prepare_passage_inputs_spatialDLPFC.R`: creates per-section
  PASSAGE input RDS files and discovers available cell-type proportion columns.
- `scripts/01_prepare_passage_inputs_spatialDLPFC.sbatch`: Biowulf batch wrapper.
- `analysis_plan_spatialDLPFC.md`: analysis plan and QC decisions.
