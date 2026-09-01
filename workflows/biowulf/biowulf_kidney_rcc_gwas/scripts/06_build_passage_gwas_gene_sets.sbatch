#!/bin/bash
#SBATCH --job-name=passage_kidney_gwasets
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=01:00:00
#SBATCH --output=logs/06_gwas_gene_sets.%j.out
#SBATCH --error=logs/06_gwas_gene_sets.%j.err

set -euo pipefail
module load R
ROOT="${1:-$PWD}"
mkdir -p "$ROOT/logs" "$ROOT/results/passage_gwas_inputs"
Rscript "$ROOT/scripts/06_build_passage_gwas_gene_sets.R" "$ROOT"
