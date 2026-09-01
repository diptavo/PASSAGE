#!/bin/bash
#SBATCH --job-name=passage_kid_pathsum
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=00:45:00
#SBATCH --output=logs/03_pathsum.%j.out
#SBATCH --error=logs/03_pathsum.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
mkdir -p "$ROOT/logs"
Rscript "$ROOT/scripts/03_aggregate_pathway_results.R" "$ROOT"
