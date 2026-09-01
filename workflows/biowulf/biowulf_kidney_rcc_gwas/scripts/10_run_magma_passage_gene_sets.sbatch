#!/bin/bash
#SBATCH --job-name=passage_kidney_magma_sets
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=04:00:00
#SBATCH --output=logs/10_magma_sets.%j.out
#SBATCH --error=logs/10_magma_sets.%j.err

set -euo pipefail
module load magma/1.10

ROOT="${1:-$PWD}"
GENE_LOC="${2:-/data/Dutta_lab/tools/NCBI38.gene.loc}"
GENE_DIR="$ROOT/results/magma_gene"
SET_DIR="$ROOT/results/passage_gwas_inputs"
ENTREZ_SET_DIR="$SET_DIR/entrez_gmt"
OUTDIR="$ROOT/results/magma_passage_sets"
mkdir -p "$ROOT/logs" "$OUTDIR"

SETS=(
  "passage_spasset_driver_sets.gmt"
  "passage_significant_scorez_driver_sets.gmt"
)

test -s "$GENE_LOC" || { echo "Missing gene location file: $GENE_LOC" >&2; exit 2; }
python3 "$ROOT/scripts/10_convert_gmt_symbols_to_entrez.py" "$GENE_LOC" "$ENTREZ_SET_DIR" \
  "$SET_DIR/${SETS[0]}" "$SET_DIR/${SETS[1]}"

for pheno in CC2 PRCC RCC; do
  GENE_RESULTS="$GENE_DIR/${pheno}.genes.raw"
  test -s "$GENE_RESULTS" || { echo "Missing $GENE_RESULTS" >&2; exit 2; }
  for set_file in "${SETS[@]}"; do
    SET_PATH="$ENTREZ_SET_DIR/${set_file%.gmt}.entrez.gmt"
    test -s "$SET_PATH" || { echo "Missing $SET_PATH" >&2; exit 2; }
    stem="${set_file%.gmt}.entrez"
    magma --gene-results "$GENE_RESULTS" \
      --set-annot "$SET_PATH" \
      --out "$OUTDIR/${pheno}.${stem}"
  done
done

find "$OUTDIR" -maxdepth 1 -type f | sort
