#!/bin/bash
#SBATCH --job-name=passage_kid_drvman
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --output=logs/04_drvman.%j.out
#SBATCH --error=logs/04_drvman.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
MAX_PATHWAYS="${2:-80}"
TOP_PER_COLLECTION_STAT="${3:-8}"
Q_CUT="${4:-0.10}"
mkdir -p "$ROOT/logs" "$ROOT/metadata"
Rscript "$ROOT/scripts/04_make_driver_manifest.R" "$ROOT" "$MAX_PATHWAYS" "$TOP_PER_COLLECTION_STAT" "$Q_CUT"
