#!/bin/bash
#SBATCH --job-name=passage_dlpfc_download
#SBATCH --cpus-per-task=4
#SBATCH --mem=48g
#SBATCH --time=12:00:00
#SBATCH --output=logs/00_download_spatialDLPFC.%j.out
#SBATCH --error=logs/00_download_spatialDLPFC.%j.err

set -euo pipefail

module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/data" "$ROOT/refs" "$ROOT/cache"
Rscript "$ROOT/scripts/00_download_spatialDLPFC.R" "$ROOT"
