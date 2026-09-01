#!/usr/bin/env python3

import collections
import csv
import math
import os
import statistics
from pathlib import Path


ROOT = Path(os.environ.get("PASSAGE_ROOT", ".")).resolve()


def read(path):
    with open(ROOT / path, newline="") as handle:
        return list(csv.DictReader(handle))


def fnum(x):
    try:
        if x in ("", "NA", "NaN", "nan", None):
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def finite(vals):
    return [v for v in vals if math.isfinite(v)]


def median(vals):
    vals = finite(vals)
    return statistics.median(vals) if vals else float("nan")


def fmt(x):
    if not math.isfinite(x):
        return "NA"
    if abs(x) < 1e-3:
        return f"{x:.2e}"
    return f"{x:.4f}"


def pname(x):
    return x.replace("HALLMARK_", "")


def main():
    h = read("results/passage_kidney_summary/all_sample_h1_h3_pathways.csv")
    c = read("results/passage_kidney_summary/all_sample_competitive_h3.csv")
    cs = read("results/passage_kidney_summary/competitive_h3_pathway_summary.csv")
    dr = read("results/passage_gwas_inputs/passage_spasset_driver_gene_long.csv")

    print("## SELF_CONTAINED_COUNTS")
    for ref in sorted({r["reference"] for r in h}):
        z = [r for r in h if r["reference"] == ref]
        print(
            ref,
            "H1_FDR05", sum(fnum(r["fdr_H1"]) <= 0.05 for r in z),
            "H3_FDR05", sum(fnum(r["fdr_H3"]) <= 0.05 for r in z),
            "rows", len(z),
            "median_reduction", fmt(median(fnum(r["celltype_adjusted_reduction"]) for r in z)),
            "median_propSV", fmt(median(fnum(r["mean_propSV"]) for r in z)),
            "median_PSVS", fmt(median(fnum(r["PSVS_range"]) for r in z)),
        )

    print("\n## COMPETITIVE_COUNTS_BY_REF_STAT")
    for ref in sorted({r["reference"] for r in c}):
        for stat in ["score_z", "score_robust_z"]:
            z = [r for r in c if r["reference"] == ref and r["statistic"] == stat]
            print(
                ref, stat,
                "raw_p05", sum(fnum(r["competitive_score_p"]) <= 0.05 for r in z),
                "global_FDR05", sum(fnum(r["competitive_score_fdr_global"]) <= 0.05 for r in z),
                "empirical_FDR05", sum(fnum(r["competitive_score_empirical_leave_sample_fdr"]) <= 0.05 for r in z),
                "median_cEPSV", fmt(median(fnum(r["cEPSV"]) for r in z)),
                "median_coh", fmt(median(fnum(r["coherence_pc1"]) for r in z)),
            )

    print("\n## TOP_SUMMARY_SCOREZ")
    for ref in sorted({r["reference"] for r in cs}):
        z = [r for r in cs if r["reference"] == ref and r["statistic"] == "score_z"]
        z = sorted(z, key=lambda r: fnum(r["fisher_p"]))[:15]
        print("\nREF", ref)
        for r in z:
            print(
                pname(r["pathway"]),
                "n_raw_p05=" + r["n_raw_p05"],
                "n_fdr05=" + r["n_raw_fdr05"],
                "n_empFDR05=" + r["n_empirical_fdr05"],
                "median_p=" + fmt(fnum(r["median_p"])),
                "fisher_p=" + fmt(fnum(r["fisher_p"])),
                "fisher_fdr=" + fmt(fnum(r["fisher_fdr"])),
                "cEPSV=" + fmt(fnum(r["median_cEPSV"])),
                "coh=" + fmt(fnum(r["median_coherence_pc1"])),
            )

    print("\n## TOP_PER_SAMPLE_SCOREZ")
    for ref in sorted({r["reference"] for r in c}):
        print("\nREF", ref)
        for sample in sorted({r["spatial_sample"] for r in c}):
            z = [
                r for r in c
                if r["reference"] == ref and r["spatial_sample"] == sample and r["statistic"] == "score_z"
            ]
            z = sorted(z, key=lambda r: fnum(r["competitive_score_p"]))[:8]
            parts = [
                f"{pname(r['pathway'])}:p={fmt(fnum(r['competitive_score_p']))}/FDR={fmt(fnum(r['competitive_score_fdr_global']))}"
                for r in z
            ]
            print(sample, "; ".join(parts))

    print("\n## DRIVERS")
    sig_keys = {
        (r["reference"], r["pathway"])
        for r in cs
        if r["statistic"] == "score_z" and fnum(r["fisher_fdr"]) <= 0.05
    }
    for ref in sorted({r["reference"] for r in dr}):
        dd = [r for r in dr if r["reference"] == ref and (r["reference"], r["pathway"]) in sig_keys]
        cnt = collections.Counter(r["gene"] for r in dd)
        print(
            ref,
            "n_driver_rows", len(dd),
            "n_unique_genes", len(cnt),
            "top_genes", ";".join([f"{g}:{n}" for g, n in cnt.most_common(20)]),
        )

    print("\n## GWAS_FILES")
    for path in [
        "refs/gwas/magma_pval_inputs/magma_pval_input_manifest.csv",
        "results/passage_gwas_inputs/README.md",
    ]:
        print("__" + path + "__")
        with open(ROOT / path) as handle:
            for i, line in enumerate(handle):
                if i > 20:
                    break
                print(line.rstrip())


if __name__ == "__main__":
    main()
