#!/bin/bash
#SBATCH --job-name=passage_kidney_gwasfmt
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_gwas_format.%j.out
#SBATCH --error=logs/05_gwas_format.%j.err

set -euo pipefail
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/refs/gwas/magma_pval_inputs"
export PASSAGE_ROOT="$ROOT"
python3 "$ROOT/scripts/05_prepare_gwas_sumstats_for_magma.py"
