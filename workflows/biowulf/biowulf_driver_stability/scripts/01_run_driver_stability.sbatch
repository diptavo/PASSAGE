#!/bin/bash
#SBATCH --job-name=passage_drv_stab
#SBATCH --cpus-per-task=1
#SBATCH --mem=28g
#SBATCH --time=00:55:00
#SBATCH --output=logs/01_drv_stab.%A_%a.out
#SBATCH --error=logs/01_drv_stab.%A_%a.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
N_BOOT="${2:-40}"
N_NULL="${3:-5}"
TOP_K="${4:-10}"
mkdir -p "$ROOT/logs" "$ROOT/results/driver_stability"
Rscript "$ROOT/scripts/01_run_driver_stability.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$N_BOOT" "$N_NULL" "$TOP_K"
