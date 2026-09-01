#!/bin/bash
#SBATCH --job-name=passage_cond_mod
#SBATCH --cpus-per-task=1
#SBATCH --mem=48g
#SBATCH --time=02:00:00
#SBATCH --output=logs/03_modules.%j.out
#SBATCH --error=logs/03_modules.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
Q_CUT="${2:-0.10}"
JACCARD="${3:-0.25}"
OVERLAP="${4:-0.50}"
MAX_CANDIDATES="${5:-2500}"
PRIMARY_STATS="${6:-score_z,score_z_robust,CSPS,GSPS,HCPS}"
mkdir -p "$ROOT/logs"
Rscript "$ROOT/scripts/03_condense_pathway_modules.R" "$ROOT" "$Q_CUT" "$JACCARD" "$OVERLAP" "$MAX_CANDIDATES" "$PRIMARY_STATS"
