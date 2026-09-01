#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path


def load_symbol_map(gene_loc):
    mapping = defaultdict(list)
    with open(gene_loc) as handle:
        for line in handle:
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                continue
            entrez, symbol = parts[0], parts[5].upper()
            if entrez not in mapping[symbol]:
                mapping[symbol].append(entrez)
    return mapping


def convert_gmt(in_gmt, out_gmt, mapping):
    stats = []
    with open(in_gmt) as fin, open(out_gmt, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, desc, genes = parts[0], parts[1], parts[2:]
            entrez = []
            matched_symbols = 0
            for gene in genes:
                ids = mapping.get(gene.upper(), [])
                if ids:
                    matched_symbols += 1
                entrez.extend(ids)
            seen = set()
            entrez_unique = []
            for gene_id in entrez:
                if gene_id not in seen:
                    seen.add(gene_id)
                    entrez_unique.append(gene_id)
            if len(entrez_unique) >= 2:
                fout.write("\t".join([name, desc] + entrez_unique) + "\n")
            stats.append({
                "set": name,
                "input_symbols": len(genes),
                "matched_symbols": matched_symbols,
                "output_entrez": len(entrez_unique),
                "kept": int(len(entrez_unique) >= 2),
            })
    return stats


def main():
    if len(sys.argv) < 4:
        raise SystemExit("usage: convert_gmt_symbols_to_entrez.py GENE_LOC OUT_DIR GMT [GMT...]")
    gene_loc = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping = load_symbol_map(gene_loc)
    all_stats = []
    for gmt in map(Path, sys.argv[3:]):
        out_gmt = out_dir / (gmt.stem + ".entrez.gmt")
        stats = convert_gmt(gmt, out_gmt, mapping)
        for row in stats:
            row["source_gmt"] = str(gmt)
            row["out_gmt"] = str(out_gmt)
        all_stats.extend(stats)
        print(out_gmt)
    with open(out_dir / "conversion_manifest.csv", "w", newline="") as fout:
        fields = ["source_gmt", "out_gmt", "set", "input_symbols", "matched_symbols", "output_entrez", "kept"]
        writer = csv.DictWriter(fout, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_stats)


if __name__ == "__main__":
    main()
