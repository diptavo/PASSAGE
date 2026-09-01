#!/bin/bash
#SBATCH --job-name=passage_dlpfc_agg
#SBATCH --cpus-per-task=2
#SBATCH --mem=24g
#SBATCH --time=2:00:00
#SBATCH --output=logs/04_aggregate.%j.out
#SBATCH --error=logs/04_aggregate.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/results/passage_dlpfc_summary"
Rscript "$ROOT/scripts/04_aggregate_passage_dlpfc.R" "$ROOT"
