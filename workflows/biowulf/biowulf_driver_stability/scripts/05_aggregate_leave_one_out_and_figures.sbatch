#!/bin/bash
#SBATCH --job-name=passage_drv_fig
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=00:25:00
#SBATCH --output=logs/05_drv_fig.%j.out
#SBATCH --error=logs/05_drv_fig.%j.err

set -euo pipefail
module load R

ROOT="${1:?root required}"
Rscript "$ROOT/scripts/05_aggregate_leave_one_out_and_figures.R" "$ROOT"
