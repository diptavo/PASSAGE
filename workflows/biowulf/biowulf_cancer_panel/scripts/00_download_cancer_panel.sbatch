#!/bin/bash
#SBATCH --job-name=passage_cancer_dl
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=08:00:00
#SBATCH --output=logs/00_download.%j.out
#SBATCH --error=logs/00_download.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/data" "$ROOT/refs"
Rscript "$ROOT/scripts/00_download_cancer_panel.R" "$ROOT"
