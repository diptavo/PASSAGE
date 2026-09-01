#!/bin/bash
#SBATCH --job-name=passage_cancer_comp
#SBATCH --cpus-per-task=2
#SBATCH --mem=80g
#SBATCH --time=18:00:00
#SBATCH --array=1-8
#SBATCH --output=logs/03_competitive_h3.%A_%a.out
#SBATCH --error=logs/03_competitive_h3.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
B="${2:-199}"
CORES="${3:-2}"
mkdir -p "$ROOT/logs" "$ROOT/results/passage_cancer_competitive_h3"
Rscript "$ROOT/scripts/03_run_passage_cancer_competitive_h3.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$B" "$CORES"
