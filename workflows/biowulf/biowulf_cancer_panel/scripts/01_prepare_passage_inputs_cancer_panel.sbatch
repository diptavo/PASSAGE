#!/bin/bash
#SBATCH --job-name=passage_cancer_prep
#SBATCH --cpus-per-task=2
#SBATCH --mem=80g
#SBATCH --time=06:00:00
#SBATCH --output=logs/01_prepare_inputs.%j.out
#SBATCH --error=logs/01_prepare_inputs.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/data/passage_inputs"
Rscript "$ROOT/scripts/01_prepare_passage_inputs_cancer_panel.R" "$ROOT"
