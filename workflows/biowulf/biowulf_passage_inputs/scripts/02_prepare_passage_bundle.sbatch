#!/usr/bin/env bash
#SBATCH --job-name=passage_inputs
#SBATCH --partition=norm
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48g
#SBATCH --gres=lscratch:50
#SBATCH --output=logs/02_prepare_passage_bundle.%A_%a.out
#SBATCH --error=logs/02_prepare_passage_bundle.%A_%a.err

set -euo pipefail
root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"
module load R
Rscript "${root}/scripts/02_prepare_passage_bundle.R" "${root}" "${SLURM_ARRAY_TASK_ID}"
