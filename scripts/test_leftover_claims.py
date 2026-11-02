#!/usr/bin/env python3
"""Answer leftover neuron/aging questions from processed GEO only."""
from __future__ import annotations

import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path("/home/noahtjones/hu-rnmp-pipeline")
PUB = ROOT / "data/public"
OUT = ROOT / "results/_cross"

GENES = {
    # interferon / cGAS–STING
    "CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "IRF7",
    "CXCL10", "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "IFIT1", "IFIT2", "IFIT3", "RSAD2", "DDX58", "IFIH1",
    "STAT1", "STAT2", "IRF9", "IFI44", "IFI27",
    # FTD / tau / paper's OPN
    "MAPT", "GRN", "C9orf72", "TARDBP", "VCP", "CHMP2B", "SPP1",
    "TMEM106B", "GSK3B", "CDK5", "MARK2", "MARK3",
    # RER
    "RNASEH2A", "RNASEH2B", "RNASEH2C", "RNASEH1", "FEN1", "POLD1",
    "LIG1", "TOP1", "TDP1",
}

ISG = [
    "CGAS", "STING1", "TBK1", "IRF3", "IRF7", "CXCL10", "ISG15", "MX1",
    "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2", "IFIT3", "RSAD2", "DDX58",
    "IFIH1", "STAT1", "STAT2", "IRF9", "IFI44", "IFI27", "MX2",
]
FTD = ["MAPT", "GRN", "C9orf72", "TARDBP", "VCP", "CHMP2B", "SPP1",
       "TMEM106B", "GSK3B", "CDK5", "MARK2"]
RER = ["RNASEH2A", "RNASEH2B", "RNASEH2C", "RNASEH1", "FEN1", "POLD1",
       "LIG1", "TOP1", "TDP1"]

# from GSE230519 series matrix Sample_title
FTD_PM = {"HO001", "GO004", "GO006", "HO008", "HO010", "HO012"}
CTRL_PM = {"HO002", "GO003", "GO005", "HO007", "HO009", "HO011", "HO013"}


def mwu(a, b):
    a = [x for x in a if x == x]
    b = [x for x in b if x == x]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def lfc(num, den):
    if num != num or den != den or num <= 0 or den <= 0:
        return float("nan")
    return math.log2(num / den)


def bed_sums(path: Path):
    """Return (total, by_chrom) weighted by interval length."""
    tot = 0.0
    by = defaultdict(float)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line or line[0] in "#t":
                continue
            p = line.split()
            if len(p) < 4:
                continue
            try:
                s, e, v = int(p[1]), int(p[2]), float(p[3])
            except ValueError:
                continue
            w = v * max(e - s, 1)
            tot += w
            by[p[0]] += w
    return tot, by


def test_nuclear_brain():
    d = PUB / "nuclear_brain"
    by = defaultdict(lambda: {"KOH": 0.0, "KCl": 0.0})
    chrom_by = defaultdict(lambda: defaultdict(lambda: {"KOH": 0.0, "KCl": 0.0}))
    pat = re.compile(r"-(?P<age>[oy]\d)-mouse-brain-(?P<chem>KCl|KOH)-")
    for p in sorted(d.glob("*.gz")):
        m = pat.search(p.name)
        if not m:
            continue
        tot, byc = bed_sums(p)
        key = m.group("age")
        chem = m.group("chem")
        by[key][chem] += tot
        for c, v in byc.items():
            chrom_by[key][c][chem] += v
    rows = []
    for aid, dlt in sorted(by.items()):
        koh, kcl = dlt["KOH"], dlt["KCl"]
        rows.append({
            "age": "old" if aid[0] == "o" else "young",
            "animal": aid,
            "koh": koh,
            "kcl": kcl,
            "koh_over_kcl": koh / kcl if kcl else float("nan"),
        })
    old = [r["koh_over_kcl"] for r in rows if r["age"] == "old"]
    yng = [r["koh_over_kcl"] for r in rows if r["age"] == "young"]
    # chrom-level old/young using pooled animals
    chrom_sum = {}
    for aid, cmap in chrom_by.items():
        age = "old" if aid[0] == "o" else "young"
        for c, dlt in cmap.items():
            chrom_sum.setdefault(c, {"old": [], "young": []})
            if dlt["KCl"] > 0:
                chrom_sum[c][age].append(dlt["KOH"] / dlt["KCl"])
    top = []
    for c, ages in chrom_sum.items():
        if len(ages["old"]) == 3 and len(ages["young"]) == 3:
            om, ym = float(np.mean(ages["old"])), float(np.mean(ages["young"]))
            top.append((c, om, ym, om / ym if ym else float("nan")))
    top.sort(key=lambda x: -abs(math.log2(x[3])) if x[3] == x[3] and x[3] > 0 else 0)
    return {
        "per_animal": rows,
        "old_mean": float(np.mean(old)),
        "young_mean": float(np.mean(yng)),
        "old_over_young": float(np.mean(old) / np.mean(yng)),
        "mwu_p": mwu(old, yng),
        "old": old,
        "young": yng,
        "chrom_most_shifted": [
            {"chrom": c, "old": om, "young": ym, "old_over_young": r}
            for c, om, ym, r in top[:12]
        ],
    }


def test_ftd_postmortem():
    """Pseudobulk (sum across cells) then CPM for predeclared genes."""
    samples = {}
    for p in sorted((PUB / "gse230519_pm").glob("*.gz")):
        m = re.search(r"_(HO\d+|GO\d+)\.txt", p.name)
        if not m:
            continue
        sid = m.group(1)
        n_cells = 0
        sums = {g: 0.0 for g in GENES}
        lib = 0.0
        with gzip.open(p, "rt") as fh:
            header = fh.readline()
            n_cells = header.count(" ") + 1  # quoted barcodes space-separated
            for line in fh:
                # "GENE" 0 1 0 ...
                if not line.startswith('"'):
                    continue
                q = line.find('"', 1)
                gene = line[1:q]
                rest = line[q + 1 :]
                # fast path: skip if gene not wanted, but still need lib size
                # lib = sum of ALL UMIs
                vals = rest.split()
                try:
                    total = sum(int(x) for x in vals)
                except ValueError:
                    continue
                lib += total
                if gene in sums:
                    sums[gene] += total
        samples[sid] = {"n_cells": n_cells, "lib": lib, "sums": sums}

    def cpm(sid, g):
        s = samples[sid]
        return 1e6 * s["sums"][g] / s["lib"] if s["lib"] else float("nan")

    def set_stats(genes):
        recs = []
        fcs = []
        for g in genes:
            ft = [cpm(s, g) for s in FTD_PM if s in samples]
            ct = [cpm(s, g) for s in CTRL_PM if s in samples]
            fc = lfc(float(np.nanmean(ft)), float(np.nanmean(ct)))
            recs.append({
                "gene": g,
                "FTD_mean_CPM": float(np.nanmean(ft)),
                "Ctrl_mean_CPM": float(np.nanmean(ct)),
                "log2fc": fc,
                "mwu_p": mwu(ft, ct),
                "n_FTD": len(ft),
                "n_Ctrl": len(ct),
            })
            if fc == fc:
                fcs.append(fc)
        return {
            "median_log2fc": float(np.median(fcs)) if fcs else float("nan"),
            "n_up_0.5": sum(1 for x in fcs if x >= 0.5),
            "n_down_0.5": sum(1 for x in fcs if x <= -0.5),
            "genes": recs,
        }

    return {
        "n_samples": {s: samples[s]["n_cells"] for s in samples},
        "libsizes": {s: samples[s]["lib"] for s in samples},
        "ISG": set_stats(ISG),
        "FTD_tau": set_stats(FTD),
        "RER": set_stats(RER),
        "SPP1": set_stats(["SPP1"]),
        "STING1": set_stats(["STING1"]),
        "MAPT": set_stats(["MAPT"]),
    }


def main():
    blob = {
        "T1b_nuclear_brain": test_nuclear_brain(),
        "T2b_GSE230519_postmortem": test_ftd_postmortem(),
        "T3b_PRJNA941970": {
            "call": "SAME_AS_GSE297046",
            "note": (
                "ENA run titles FS303/FS307/FS310/FS326 match GSE297046 "
                "ribome bedgraphs. No separate RNASEH2A-KO mt table. "
                "FASTQ-only besides 350 MB bigWigs and 500 kb bins already scored."
            ),
        },
        "T4b_PXD069495": {
            "call": "UNTESTABLE_NO_SMALL_TABLE",
            "note": (
                "PRIDE PARTIAL deposit: ~12 GB raw zip per sample; FDR xlsx "
                "25–57 MB each (per-sample search output). No project-level "
                "PTM site table. Tau PTM ≠ RNase H2 / rNMP."
            ),
        },
    }
    t1 = blob["T1b_nuclear_brain"]
    t2 = blob["T2b_GSE230519_postmortem"]
    blob["verdicts"] = {
        "T1b_nuclear_brain_age": {
            "old_over_young": t1["old_over_young"],
            "p": t1["mwu_p"],
            "call": (
                "NOT_SUPPORTED_AT_N3"
                if t1["old_over_young"] <= 1.05
                else "WEAK_OR_UNDERPOWERED"
                if t1["old_over_young"] < 1.3
                else "SUPPORTED_DIRECTION"
            ),
        },
        "T2b_FTD_brain_ISG_OPN": {
            "ISG_median_log2fc": t2["ISG"]["median_log2fc"],
            "SPP1_log2fc": t2["SPP1"]["genes"][0]["log2fc"],
            "STING1_log2fc": t2["STING1"]["genes"][0]["log2fc"],
            "MAPT_log2fc": t2["MAPT"]["genes"][0]["log2fc"],
            "call": "see per-gene; OPN is the paper claim",
        },
    }
    (OUT / "TEST_RESULTS_LEFTOVER.json").write_text(json.dumps(blob, indent=2, default=str))
    print(json.dumps(blob["verdicts"], indent=2, default=str))
    print("nuclear animals", t1["per_animal"])
    print("SPP1", t2["SPP1"])
    print("STING1", t2["STING1"])
    print("MAPT", t2["MAPT"])
    print("ISG genes")
    for g in t2["ISG"]["genes"]:
        print(g)
    print("FTD genes")
    for g in t2["FTD_tau"]["genes"]:
        print(g)
    print("RER genes")
    for g in t2["RER"]["genes"]:
        print(g)


if __name__ == "__main__":
    main()
