#!/usr/bin/env bash
#SBATCH --job-name=passage_validate_inputs
#SBATCH --partition=norm
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48g
#SBATCH --gres=lscratch:20
#SBATCH --output=logs/03_validate_inputs.%j.out
#SBATCH --error=logs/03_validate_inputs.%j.err

set -euo pipefail
root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"
module load R
Rscript "${root}/scripts/03_validate_inputs.R" "${root}"
