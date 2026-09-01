#!/bin/bash
#SBATCH --job-name=passage_cond_mtest
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=03:00:00
#SBATCH --output=logs/05_module.%A_%a.out
#SBATCH --error=logs/05_module.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
B="${2:-999}"
mkdir -p "$ROOT/logs" "$ROOT/results/module_testing"
Rscript "$ROOT/scripts/05_run_module_task.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$B"
