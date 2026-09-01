#!/bin/bash

set -euo pipefail

ROOT="${1:?Usage: submit_pipeline.sh ROOT PASSAGE_CODE [N_PERM] [B]}"
PASSAGE_CODE="${2:?Provide the PASSAGE repository path on Biowulf}"
N_PERM="${3:-999}"
B="${4:-9999}"

mkdir -p "$ROOT/logs" "$ROOT/r_libs"
cd "$ROOT"
export PASSAGE_CODE R_LIBS_USER="$ROOT/r_libs"

download_job=$(sbatch --parsable scripts/00_download_kidney_rcc_spatial_refs.sbatch "$ROOT")
msigdb_job=$(sbatch --parsable scripts/01_prepare_hallmark_cache.sbatch "$ROOT")
prepare_job=$(sbatch --parsable \
  --dependency="afterok:${download_job}:${msigdb_job}" \
  scripts/01_prepare_passage_inputs_kidney_rcc.sbatch "$ROOT")
h1_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" \
  --array=1-6 scripts/02_run_passage_kidney_h1_h3.sbatch \
  "$ROOT" "$N_PERM" 4)
competitive_job=$(sbatch --parsable --dependency="afterok:${h1_job}" \
  --array=1-6 scripts/03_run_passage_kidney_competitive_h3.sbatch \
  "$ROOT" "$B" 2)
aggregate_job=$(sbatch --parsable \
  --dependency="afterok:${h1_job}:${competitive_job}" \
  scripts/04_aggregate_passage_kidney.sbatch "$ROOT")

printf 'download=%s\nmsigdb=%s\nprepare=%s\nh1_h3=%s\ncompetitive_h3=%s\naggregate=%s\n' \
  "$download_job" "$msigdb_job" "$prepare_job" "$h1_job" \
  "$competitive_job" "$aggregate_job"
