#!/usr/bin/env python3

import csv
import os
import re
from pathlib import Path


ROOT = Path(os.environ.get("PASSAGE_ROOT", ".")).resolve()
OUTDIR = ROOT / "results" / "magma_passage_sets"
SUMMARY = ROOT / "results" / "magma_passage_sets_summary"
SUMMARY.mkdir(parents=True, exist_ok=True)


def read_magma_table(path):
    rows = []
    with open(path) as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if header is None:
                header = parts
                continue
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
    return rows


all_rows = []
for path in sorted(OUTDIR.glob("*.gsa.out")):
    name = path.name.replace(".gsa.out", "")
    phenotype, set_name = name.split(".", 1)
    for row in read_magma_table(path):
        row["phenotype"] = phenotype
        row["set_collection"] = set_name
        all_rows.append(row)

if all_rows:
    columns = ["phenotype", "set_collection"] + [c for c in all_rows[0].keys() if c not in ("phenotype", "set_collection")]
else:
    columns = ["phenotype", "set_collection"]

out_csv = SUMMARY / "magma_passage_gene_set_results.csv"
with open(out_csv, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

print(out_csv)
