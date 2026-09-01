#!/bin/bash
#SBATCH --job-name=passage_drv_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=00:20:00
#SBATCH --output=logs/02_drv_sum.%j.out
#SBATCH --error=logs/02_drv_sum.%j.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
mkdir -p "$ROOT/logs" "$ROOT/results/driver_stability_summary"
Rscript "$ROOT/scripts/02_aggregate_driver_stability.R" "$ROOT"
