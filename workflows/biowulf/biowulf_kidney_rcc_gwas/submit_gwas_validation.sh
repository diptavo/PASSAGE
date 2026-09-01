#!/bin/bash

set -euo pipefail

ROOT="${1:?Usage: submit_gwas_validation.sh ROOT [BFILE] [GENE_LOC]}"
BFILE="${2:-/data/Dutta_lab/REF/EUR}"
GENE_LOC="${3:-/data/Dutta_lab/tools/NCBI38.gene.loc}"
ANNOT="$ROOT/refs/gwas/magma_annotation/EUR_NCBI38.genes.annot"

mkdir -p "$ROOT/logs"
cd "$ROOT"

format_job=$(sbatch --parsable scripts/05_prepare_gwas_sumstats_for_magma.sbatch "$ROOT")
sets_job=$(sbatch --parsable scripts/06_build_passage_gwas_gene_sets.sbatch "$ROOT")
annot_job=$(sbatch --parsable scripts/08_prepare_magma_annotation.sbatch \
  "$ROOT" "$BFILE" "$GENE_LOC")
gene_job=$(sbatch --parsable \
  --dependency="afterok:${format_job}:${annot_job}" --array=1-3 \
  scripts/09_run_magma_gene_analysis.sbatch "$ROOT" "$BFILE" "$ANNOT")
test_job=$(sbatch --parsable \
  --dependency="afterok:${sets_job}:${gene_job}" \
  scripts/10_run_magma_passage_gene_sets.sbatch "$ROOT" "$GENE_LOC")
summary_job=$(sbatch --parsable --dependency="afterok:${test_job}" \
  scripts/11_summarize_magma_passage_sets.sbatch "$ROOT")

printf 'format_sumstats=%s\nbuild_sets=%s\nannotation=%s\ngene_analysis=%s\nset_tests=%s\nsummary=%s\n' \
  "$format_job" "$sets_job" "$annot_job" "$gene_job" "$test_job" \
  "$summary_job"
