#!/bin/bash
#SBATCH --job-name=passage_kid_driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=36g
#SBATCH --time=06:00:00
#SBATCH --output=logs/05_driver.%A_%a.out
#SBATCH --error=logs/05_driver.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
N_BOOT="${2:-30}"
N_NULL="${3:-9}"
TOP_K="${4:-10}"
mkdir -p "$ROOT/logs" "$ROOT/results/driver_stability"
Rscript "$ROOT/scripts/05_run_driver_stability.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$N_BOOT" "$N_NULL" "$TOP_K"
