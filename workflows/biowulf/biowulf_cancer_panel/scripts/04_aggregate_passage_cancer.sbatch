#!/bin/bash
#SBATCH --job-name=passage_cancer_agg
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --output=logs/04_aggregate.%j.out
#SBATCH --error=logs/04_aggregate.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/results/passage_cancer_summary"
Rscript "$ROOT/scripts/04_aggregate_passage_cancer.R" "$ROOT"
