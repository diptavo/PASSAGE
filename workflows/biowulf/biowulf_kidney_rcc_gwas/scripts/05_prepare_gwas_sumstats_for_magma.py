#!/usr/bin/env python3

import csv
import gzip
import os
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path, "rt")


def main():
    root = Path(os.environ.get("PASSAGE_ROOT", ".")).resolve()
    manifest = root / "refs" / "gwas" / "rcc_gwas_sumstats_manifest.csv"
    out_dir = root / "refs" / "gwas" / "magma_pval_inputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(manifest, newline="")))
    status = []
    for row in rows:
        phenotype = row["phenotype"]
        src = Path(row["file"])
        out = out_dir / f"{phenotype}.magma.pval.tsv"
        n = 0
        kept = 0
        with open_text(src) as fin, open(out, "w", newline="") as fout:
            reader = csv.DictReader(fin, delimiter="\t")
            writer = csv.writer(fout, delimiter="\t")
            writer.writerow(["SNP", "P", "N", "CHR", "BP", "A1", "A2", "BETA", "SE"])
            for rec in reader:
                n += 1
                snp = rec.get("rsid") or rec.get("variant_id")
                p = rec.get("p_value")
                ns = rec.get("n_samples")
                if not snp or not p or p in ("NA", "nan", ""):
                    continue
                writer.writerow([
                    snp, p, ns,
                    rec.get("chromosome", ""),
                    rec.get("base_pair_location", ""),
                    rec.get("effect_allele", ""),
                    rec.get("other_allele", ""),
                    rec.get("beta", ""),
                    rec.get("standard_error", ""),
                ])
                kept += 1
        status.append({
            "phenotype": phenotype,
            "source": str(src),
            "out": str(out),
            "input_rows": n,
            "kept_rows": kept,
        })
    with open(out_dir / "magma_pval_input_manifest.csv", "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=["phenotype", "source", "out", "input_rows", "kept_rows"])
        writer.writeheader()
        writer.writerows(status)


if __name__ == "__main__":
    main()
