#!/bin/bash
#SBATCH --job-name=passage_kid_manifest
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --output=logs/01_manifest.%j.out
#SBATCH --error=logs/01_manifest.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
CHUNK_SIZE="${2:-300}"
KIDNEY_ROOT="${3:-/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803}"
mkdir -p "$ROOT/logs" "$ROOT/metadata"
Rscript "$ROOT/scripts/01_make_pathway_manifest.R" "$ROOT" "$CHUNK_SIZE" "$KIDNEY_ROOT"
