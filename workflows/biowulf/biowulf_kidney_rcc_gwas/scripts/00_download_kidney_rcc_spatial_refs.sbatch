#!/bin/bash
#SBATCH --job-name=passage_kidney_dl
#SBATCH --cpus-per-task=2
#SBATCH --mem=24g
#SBATCH --time=08:00:00
#SBATCH --output=logs/00_download.%j.out
#SBATCH --error=logs/00_download.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/data" "$ROOT/refs"
Rscript "$ROOT/scripts/00_download_kidney_rcc_spatial_refs.R" "$ROOT"
