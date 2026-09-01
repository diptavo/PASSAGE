#!/usr/bin/env bash
#SBATCH --job-name=passage_stage_refs
#SBATCH --partition=norm
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --gres=lscratch:20
#SBATCH --output=logs/00_stage_references.%j.out
#SBATCH --error=logs/00_stage_references.%j.err

set -euo pipefail
root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"
bash "${root}/scripts/00_stage_references.sh" "${root}"
