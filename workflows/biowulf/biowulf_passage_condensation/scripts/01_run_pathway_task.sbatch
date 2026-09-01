#!/bin/bash
#SBATCH --job-name=passage_cond_path
#SBATCH --cpus-per-task=1
#SBATCH --mem=36g
#SBATCH --time=05:00:00
#SBATCH --output=logs/01_pathway.%A_%a.out
#SBATCH --error=logs/01_pathway.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
B="${2:-199}"
MIN_GENES="${3:-15}"
mkdir -p "$ROOT/logs" "$ROOT/results/pathway_testing"
Rscript "$ROOT/scripts/01_run_pathway_task.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$B" "$MIN_GENES"
