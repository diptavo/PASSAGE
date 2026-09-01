#!/bin/bash
#SBATCH --job-name=passage_path_stats
#SBATCH --cpus-per-task=1
#SBATCH --mem=28g
#SBATCH --time=00:55:00
#SBATCH --output=logs/01_path_stats.%A_%a.out
#SBATCH --error=logs/01_path_stats.%A_%a.err

set -euo pipefail
module load R

BENCH_ROOT="${1:?bench root required}"
COHORT="${2:?cohort required: breast|kidney|dlpfc}"
B="${3:-199}"
MIN_GENES="${4:-15}"
MAX_GENES="${5:-250}"

mkdir -p "$BENCH_ROOT/logs" "$BENCH_ROOT/results/pathway_new_stats/$COHORT"
Rscript "$BENCH_ROOT/scripts/01_run_pathway_new_stats.R" "$BENCH_ROOT" "$COHORT" "$SLURM_ARRAY_TASK_ID" "$B" "$MIN_GENES" "$MAX_GENES"
