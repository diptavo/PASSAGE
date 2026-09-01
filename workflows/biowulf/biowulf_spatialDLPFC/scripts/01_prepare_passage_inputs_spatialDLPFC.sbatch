#!/bin/bash
#SBATCH --job-name=passage_dlpfc_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=64g
#SBATCH --time=8:00:00
#SBATCH --output=logs/01_prepare_passage_inputs.%j.out
#SBATCH --error=logs/01_prepare_passage_inputs.%j.err

set -euo pipefail

module load R
ROOT="${1:-$PWD}"
MAX_SAMPLES="${2:-NA}"
mkdir -p "$ROOT/logs" "$ROOT/data/passage_inputs"
Rscript "$ROOT/scripts/01_prepare_passage_inputs_spatialDLPFC.R" "$ROOT" "$MAX_SAMPLES"
