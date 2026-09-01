#!/bin/bash
#SBATCH --job-name=passage_cond_msum
#SBATCH --cpus-per-task=1
#SBATCH --mem=24g
#SBATCH --time=00:45:00
#SBATCH --output=logs/06_modulesum.%j.out
#SBATCH --error=logs/06_modulesum.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
mkdir -p "$ROOT/logs"
Rscript "$ROOT/scripts/06_aggregate_module_results.R" "$ROOT"
