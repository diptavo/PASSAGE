#!/bin/bash
#SBATCH --job-name=passage_kidney_magma_gene
#SBATCH --cpus-per-task=2
#SBATCH --mem=24g
#SBATCH --time=08:00:00
#SBATCH --array=1-3
#SBATCH --output=logs/09_magma_gene.%A_%a.out
#SBATCH --error=logs/09_magma_gene.%A_%a.err

set -euo pipefail
module load magma/1.10

ROOT="${1:-$PWD}"
BFILE="${2:-/data/Dutta_lab/REF/EUR}"
ANNOT="${3:-$ROOT/refs/gwas/magma_annotation/EUR_NCBI38.genes.annot}"
MANIFEST="$ROOT/refs/gwas/magma_pval_inputs/magma_pval_input_manifest.csv"
OUTDIR="$ROOT/results/magma_gene"
mkdir -p "$ROOT/logs" "$OUTDIR"

test -s "$MANIFEST" || { echo "Missing $MANIFEST" >&2; exit 2; }
test -s "$ANNOT" || { echo "Missing $ANNOT" >&2; exit 2; }
for suffix in bed bim fam; do
  test -s "${BFILE}.${suffix}" || { echo "Missing ${BFILE}.${suffix}" >&2; exit 2; }
done

ROW=$(awk -F, -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print}' "$MANIFEST")
test -n "$ROW" || { echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
PHENO=$(printf '%s\n' "$ROW" | awk -F, '{print $1}')
PVAL=$(printf '%s\n' "$ROW" | awk -F, '{print $3}')
test -s "$PVAL" || { echo "Missing p-value input $PVAL" >&2; exit 2; }

OUT="$OUTDIR/$PHENO"
magma --bfile "$BFILE" \
  --pval "$PVAL" snp-id=SNP pval=P ncol=N \
  --gene-annot "$ANNOT" \
  --out "$OUT"

ls -lh "${OUT}.genes.raw" "${OUT}.genes.out"
