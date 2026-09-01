#!/bin/bash
#SBATCH --job-name=passage_drv_loo
#SBATCH --cpus-per-task=1
#SBATCH --mem=28g
#SBATCH --time=00:55:00
#SBATCH --output=logs/04_drv_loo.%A_%a.out
#SBATCH --error=logs/04_drv_loo.%A_%a.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
TOP_K="${2:-10}"
N_CONTROLS="${3:-10}"
mkdir -p "$ROOT/logs" "$ROOT/results/leave_one_out"
Rscript "$ROOT/scripts/04_run_leave_one_out.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$TOP_K" "$N_CONTROLS"
