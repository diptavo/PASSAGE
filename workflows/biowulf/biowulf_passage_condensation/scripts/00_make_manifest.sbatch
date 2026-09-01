#!/bin/bash
#SBATCH --job-name=passage_cond_man
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --output=logs/00_manifest.%j.out
#SBATCH --error=logs/00_manifest.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
COHORTS="${2:-kidney,breast,dlpfc}"
CHUNK_SIZE="${3:-300}"
mkdir -p "$ROOT/logs" "$ROOT/metadata"
Rscript "$ROOT/scripts/00_make_manifest.R" "$ROOT" "$COHORTS" "$CHUNK_SIZE"
