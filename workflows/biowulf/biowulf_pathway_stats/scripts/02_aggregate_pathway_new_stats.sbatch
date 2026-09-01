#!/bin/bash
#SBATCH --job-name=passage_path_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --output=logs/02_path_sum.%j.out
#SBATCH --error=logs/02_path_sum.%j.err

set -euo pipefail
module load R

BENCH_ROOT="${1:?bench root required}"
mkdir -p "$BENCH_ROOT/logs" "$BENCH_ROOT/results/pathway_new_stats_summary"
Rscript "$BENCH_ROOT/scripts/02_aggregate_pathway_new_stats.R" "$BENCH_ROOT"
