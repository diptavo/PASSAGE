#!/bin/bash
#SBATCH --job-name=passage_kid_drvsum
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH --time=00:45:00
#SBATCH --output=logs/06_drvsum.%j.out
#SBATCH --error=logs/06_drvsum.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
mkdir -p "$ROOT/logs"
Rscript "$ROOT/scripts/06_aggregate_driver_stability.R" "$ROOT"
