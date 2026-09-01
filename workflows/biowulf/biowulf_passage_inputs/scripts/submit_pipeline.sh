#!/usr/bin/env bash
set -euo pipefail

root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"
cd "${root}"
mkdir -p metadata logs

stage_job=$(sbatch --parsable scripts/00_stage_references.sbatch "${root}")
signature_job=$(sbatch --parsable --dependency="afterok:${stage_job}" scripts/01_build_reference_signatures.sbatch "${root}")
bundle_job=$(sbatch --parsable --dependency="afterok:${signature_job}" --array=1-4%4 scripts/02_prepare_passage_bundle.sbatch "${root}")
validation_job=$(sbatch --parsable --dependency="afterok:${bundle_job}" scripts/03_validate_inputs.sbatch "${root}")

cat > metadata/submitted_jobs.tsv <<EOF
stage\t${stage_job}
reference_signatures\t${signature_job}
sample_bundles\t${bundle_job}
validation\t${validation_job}
EOF
printf '%s\n' "stage=${stage_job}" "reference_signatures=${signature_job}" "sample_bundles=${bundle_job}" "validation=${validation_job}"
