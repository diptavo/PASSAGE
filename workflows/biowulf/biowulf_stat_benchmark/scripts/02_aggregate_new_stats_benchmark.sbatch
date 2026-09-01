#!/bin/bash
#SBATCH --job-name=passage_new_stats_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --output=logs/02_new_stats_summary.%j.out
#SBATCH --error=logs/02_new_stats_summary.%j.err

set -euo pipefail
module load R

BENCH_ROOT="${1:-/data/Dutta_lab/SPATH/PASSAGE_stat_benchmark_20260814}"
mkdir -p "$BENCH_ROOT/logs" "$BENCH_ROOT/results/new_stats_benchmark_summary"
Rscript "$BENCH_ROOT/scripts/02_aggregate_new_stats_benchmark.R" "$BENCH_ROOT"
