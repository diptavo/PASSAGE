#!/bin/bash
#SBATCH --job-name=passage_drv_cmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:15:00
#SBATCH --output=logs/03_drv_cmp.%j.out
#SBATCH --error=logs/03_drv_cmp.%j.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
TOP_K="${2:-10}"
Rscript "$ROOT/scripts/03_compare_driver_statistics.R" "$ROOT" "$TOP_K"
