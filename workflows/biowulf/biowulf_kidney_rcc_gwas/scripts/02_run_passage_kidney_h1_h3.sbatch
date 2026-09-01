#!/bin/bash
#SBATCH --job-name=passage_kidney_h1h3
#SBATCH --cpus-per-task=4
#SBATCH --mem=80g
#SBATCH --time=18:00:00
#SBATCH --array=1-6
#SBATCH --output=logs/02_h1_h3.%A_%a.out
#SBATCH --error=logs/02_h1_h3.%A_%a.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
N_PERM="${2:-199}"
CORES="${3:-4}"
mkdir -p "$ROOT/logs" "$ROOT/results/passage_kidney_h1_h3"
Rscript "$ROOT/scripts/02_run_passage_kidney_h1_h3.R" "$ROOT" "$SLURM_ARRAY_TASK_ID" "$N_PERM" "$CORES"
