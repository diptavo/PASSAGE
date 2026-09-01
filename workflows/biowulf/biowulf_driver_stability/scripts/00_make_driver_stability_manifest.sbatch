#!/bin/bash
#SBATCH --job-name=passage_drv_manifest
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:10:00
#SBATCH --output=logs/00_drv_manifest.%j.out
#SBATCH --error=logs/00_drv_manifest.%j.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
PATHWAY_ROOT="${2:-/data/DCEG_Dutta/PASSAGE_pathway_new_stats_20260814}"
MAX_PATHWAYS="${3:-12}"
mkdir -p "$ROOT/logs" "$ROOT/metadata"
Rscript "$ROOT/scripts/00_make_driver_stability_manifest.R" "$ROOT" "$PATHWAY_ROOT" "$MAX_PATHWAYS"
