#!/bin/bash
#SBATCH --job-name=passage_kid_gwas
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=04:00:00
#SBATCH --output=logs/07_gwas.%j.out
#SBATCH --error=logs/07_gwas.%j.err

set -euo pipefail
module load R
ROOT="${1:?root required}"
B="${2:-999}"
TOP_K="${3:-10}"
GENE_LOC="${4:-/data/Dutta_lab/tools/NCBI38.gene.loc}"
MAGMA_GENE_ROOT="${5:-/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803/results/magma_gene}"
mkdir -p "$ROOT/logs" "$ROOT/results/gwas_driver_validation"
Rscript "$ROOT/scripts/07_run_gwas_driver_enrichment.R" "$ROOT" "$B" "$TOP_K" "$GENE_LOC" "$MAGMA_GENE_ROOT"
