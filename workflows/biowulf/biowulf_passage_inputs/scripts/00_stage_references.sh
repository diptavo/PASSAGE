#!/usr/bin/env bash
set -euo pipefail

root="${1:-/data/DCEG_Dutta/PASSAGE_production_inputs_20260825}"

breast_archive="/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803/refs/raw_scRNA/breast/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
msigdb_root="/data/DCEG_Dutta/PASSAGE_condensation_multidataset_fast_20260817/refs/msigdb"
kidney_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224630/suppl"

mkdir -p "${root}/refs/breast/GSE176078" "${root}/refs/kidney/GSE224630" \
  "${root}/pathways" "${root}/metadata" "${root}/inputs" "${root}/logs"

if [[ ! -s "${root}/refs/breast/GSE176078/Wu_etal_2021_BRCA_scRNASeq/metadata.csv" ]]; then
  tar -xzf "${breast_archive}" -C "${root}/refs/breast/GSE176078"
fi

download_if_missing() {
  local url="$1"
  local destination="$2"
  if [[ ! -s "${destination}" ]]; then
    curl --fail --location --retry 5 --retry-delay 10 --output "${destination}.part" "${url}"
    mv "${destination}.part" "${destination}"
  fi
}

download_if_missing \
  "${kidney_url}/GSE224630_normalized_integrated.data.tsv.gz" \
  "${root}/refs/kidney/GSE224630/GSE224630_normalized_integrated.data.tsv.gz"
download_if_missing \
  "${kidney_url}/GSE224630_overall_metadata.tsv.gz" \
  "${root}/refs/kidney/GSE224630/GSE224630_overall_metadata.tsv.gz"

cp "${msigdb_root}/msigdb_human_pathways_filtered.rds" \
  "${root}/pathways/msigdb_human_pathways_filtered.rds"
cp "${msigdb_root}/msigdb_human_pathways_filtered_metadata.csv" \
  "${root}/pathways/msigdb_human_pathways_filtered_metadata.csv"

cat > "${root}/refs/reference_sources.tsv" <<EOF
cohort\treference_id\trole\tsource\tlocal_path
breast\tGSE176078\ttumor and breast-tissue single-cell reference\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078\t${root}/refs/breast/GSE176078
kidney\tGSE224630\tccRCC and normal-kidney single-cell reference\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630\t${root}/refs/kidney/GSE224630
pathways\tMSigDB_filtered\thuman pathway definitions\thttps://www.gsea-msigdb.org/gsea/msigdb\t${root}/pathways/msigdb_human_pathways_filtered.rds
EOF

date -u +"%Y-%m-%dT%H:%M:%SZ" > "${root}/metadata/references_staged.complete"
