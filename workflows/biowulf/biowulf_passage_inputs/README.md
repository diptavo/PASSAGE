# PASSAGE kidney and breast production inputs

This workflow prepares and validates PASSAGE inputs 1-6 for three kidney cancer Visium sections and one breast cancer Visium FFPE section.

The derived bundle root is `/data/DCEG_Dutta/PASSAGE_production_inputs_20260825`.

Stages:

1. `00_stage_references`: stage GSE176078 breast and GSE224630 kidney single-cell references plus the filtered MSigDB cache.
2. `01_build_reference_signatures`: construct donor-balanced, broad-cell-type reference signatures.
3. `02_prepare_passage_bundle`: create four self-contained sample bundles with expression, coordinates, metadata, estimated cell fractions, and adjustment designs.
4. `03_validate_inputs`: verify dimensions, identifiers, finite values, simplex constraints, design rank, pathway coverage, and deconvolution reconstruction.

All computational stages run as Biowulf batch jobs.
