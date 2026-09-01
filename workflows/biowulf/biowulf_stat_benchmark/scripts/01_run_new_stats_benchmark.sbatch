#!/bin/bash
#SBATCH --job-name=passage_new_stats
#SBATCH --cpus-per-task=1
#SBATCH --mem=28g
#SBATCH --time=00:55:00
#SBATCH --output=logs/01_new_stats.%A_%a.out
#SBATCH --error=logs/01_new_stats.%A_%a.err

set -euo pipefail
module load R

BENCH_ROOT="${1:-/data/Dutta_lab/SPATH/PASSAGE_stat_benchmark_20260814}"
COHORT="${2:?cohort required: breast|kidney|dlpfc}"
N_TARGETS="${3:-24}"
B="${4:-39}"

mkdir -p "$BENCH_ROOT/logs" "$BENCH_ROOT/results/new_stats_benchmark/$COHORT"
Rscript "$BENCH_ROOT/scripts/01_run_new_stats_benchmark.R" "$BENCH_ROOT" "$COHORT" "$SLURM_ARRAY_TASK_ID" "$N_TARGETS" "$B"
