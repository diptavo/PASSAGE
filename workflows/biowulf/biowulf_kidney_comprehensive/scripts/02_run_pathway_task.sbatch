#!/bin/bash
#SBATCH --job-name=passage_kid_path
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=04:00:00
#SBATCH --output=logs/02_pathway.%A_%a.out
#SBATCH --error=logs/02_pathway.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
B="${2:-199}"
mkdir -p "$ROOT/logs" "$ROOT/results/pathway_testing"
Rscript "$ROOT/scripts/02_run_pathway_task.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$B"
