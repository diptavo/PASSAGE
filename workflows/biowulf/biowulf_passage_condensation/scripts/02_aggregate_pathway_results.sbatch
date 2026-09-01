#!/bin/bash
#SBATCH --job-name=passage_cond_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=01:00:00
#SBATCH --output=logs/02_summary.%j.out
#SBATCH --error=logs/02_summary.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
mkdir -p "$ROOT/logs"
Rscript "$ROOT/scripts/02_aggregate_pathway_results.R" "$ROOT"
