#!/bin/bash
#SBATCH --job-name=passage_kidney_magma_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:30:00
#SBATCH --output=logs/11_magma_summary.%j.out
#SBATCH --error=logs/11_magma_summary.%j.err

set -euo pipefail
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/results/magma_passage_sets_summary"
export PASSAGE_ROOT="$ROOT"
python3 "$ROOT/scripts/11_summarize_magma_passage_sets.py"
