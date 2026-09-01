#!/bin/bash

set -euo pipefail

ROOT="${1:?Usage: submit_pipeline.sh ROOT PASSAGE_CODE [CANCERS] [N_PERM] [B]}"
PASSAGE_CODE="${2:?Provide the PASSAGE repository path on Biowulf}"
CANCERS="${3:-breast}"
N_PERM="${4:-999}"
B="${5:-9999}"

IFS=',' read -r -a CANCER_ARRAY <<< "$CANCERS"
SEEN_CANCERS=","
for cancer in "${CANCER_ARRAY[@]}"; do
  case "$cancer" in
    breast|cervical|prostate|lung) ;;
    *)
      echo "Unknown cancer '$cancer'; use breast, cervical, prostate, or lung." >&2
      exit 2
      ;;
  esac
  if [[ "$SEEN_CANCERS" == *",${cancer},"* ]]; then
    echo "Cancer '$cancer' was specified more than once." >&2
    exit 2
  fi
  SEEN_CANCERS+="${cancer},"
done
N_TASKS=$((2 * ${#CANCER_ARRAY[@]}))
if (( N_TASKS < 2 || N_TASKS > 8 )); then
  echo "CANCERS must contain one to four comma-separated workflow cancers." >&2
  exit 2
fi

mkdir -p "$ROOT/logs" "$ROOT/r_libs"
cd "$ROOT"
export PASSAGE_CODE PASSAGE_CANCERS="$CANCERS" R_LIBS_USER="$ROOT/r_libs"

download_job=$(sbatch --parsable scripts/00_download_cancer_panel.sbatch "$ROOT")
msigdb_job=$(sbatch --parsable scripts/01_prepare_hallmark_cache.sbatch "$ROOT")
prepare_job=$(sbatch --parsable \
  --dependency="afterok:${download_job}:${msigdb_job}" \
  scripts/01_prepare_passage_inputs_cancer_panel.sbatch "$ROOT")
h1_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" \
  --array="1-${N_TASKS}" scripts/02_run_passage_cancer_h1_h3.sbatch \
  "$ROOT" "$N_PERM" 4)
competitive_job=$(sbatch --parsable --dependency="afterok:${h1_job}" \
  --array="1-${N_TASKS}" scripts/03_run_passage_cancer_competitive_h3.sbatch \
  "$ROOT" "$B" 2)
aggregate_job=$(sbatch --parsable \
  --dependency="afterok:${h1_job}:${competitive_job}" \
  scripts/04_aggregate_passage_cancer.sbatch "$ROOT")

printf 'download=%s\nmsigdb=%s\nprepare=%s\nh1_h3=%s\ncompetitive_h3=%s\naggregate=%s\n' \
  "$download_job" "$msigdb_job" "$prepare_job" "$h1_job" \
  "$competitive_job" "$aggregate_job"
