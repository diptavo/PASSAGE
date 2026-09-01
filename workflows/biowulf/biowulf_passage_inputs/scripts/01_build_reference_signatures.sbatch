#!/usr/bin/env bash
#SBATCH --job-name=passage_ref_signatures
#SBATCH --partition=norm
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --output=logs/01_build_reference_signatures.%j.out
#SBATCH --error=logs/01_build_reference_signatures.%j.err

set -euo pipefail
root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"
module load R
Rscript "${root}/scripts/01_build_reference_signatures.R" "${root}"
