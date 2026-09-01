#!/bin/bash
#SBATCH --job-name=passage_cond_modman
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=00:20:00
#SBATCH --output=logs/04_modman.%j.out
#SBATCH --error=logs/04_modman.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
MAX_MODULES="${2:-250}"
mkdir -p "$ROOT/logs" "$ROOT/metadata"
Rscript "$ROOT/scripts/04_make_module_manifest.R" "$ROOT" "$MAX_MODULES"
