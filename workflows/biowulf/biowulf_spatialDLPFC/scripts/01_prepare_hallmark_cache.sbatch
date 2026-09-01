#!/bin/bash
#SBATCH --job-name=passage_dlpfc_msig
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --output=logs/01_msig.%j.out
#SBATCH --error=logs/01_msig.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/refs/msigdb"
Rscript "$ROOT/scripts/01_prepare_hallmark_cache.R" "$ROOT"
