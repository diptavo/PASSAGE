#!/bin/bash
#SBATCH --job-name=passage_kidney_magma_annot
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --output=logs/08_magma_annot.%j.out
#SBATCH --error=logs/08_magma_annot.%j.err

set -euo pipefail
module load magma/1.10

ROOT="${1:-$PWD}"
BFILE="${2:-/data/Dutta_lab/REF/EUR}"
GENE_LOC="${3:-/data/Dutta_lab/tools/NCBI38.gene.loc}"
OUTDIR="$ROOT/refs/gwas/magma_annotation"
mkdir -p "$ROOT/logs" "$OUTDIR"

for suffix in bed bim fam; do
  test -s "${BFILE}.${suffix}" || { echo "Missing ${BFILE}.${suffix}" >&2; exit 2; }
done
test -s "$GENE_LOC" || { echo "Missing gene location file: $GENE_LOC" >&2; exit 2; }

SNP_LOC="$OUTDIR/EUR.snps.loc"
ANNOT_PREFIX="$OUTDIR/EUR_NCBI38"

if [[ ! -s "$SNP_LOC" ]]; then
  awk 'BEGIN{OFS="\t"} {print $2,$1,$4}' "${BFILE}.bim" > "$SNP_LOC"
fi

if [[ ! -s "${ANNOT_PREFIX}.genes.annot" ]]; then
  magma --annotate \
    --snp-loc "$SNP_LOC" \
    --gene-loc "$GENE_LOC" \
    --out "$ANNOT_PREFIX"
fi

ls -lh "$SNP_LOC" "${ANNOT_PREFIX}.genes.annot"
