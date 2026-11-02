#!/usr/bin/env python3
"""Falsify the five neuron/aging bridges with public processed tables. No FASTQ."""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path("/home/noahtjones/hu-rnmp-pipeline")
PUB = ROOT / "data/public"
OUT = ROOT / "results/_cross"
OUT.mkdir(parents=True, exist_ok=True)

ISG_MOUSE = [
    "Cgas", "Mb21d1", "Sting1", "Tmem173", "Tbk1", "Irf3", "Irf7",
    "Cxcl10", "Isg15", "Mx1", "Oas1a", "Oas1b", "Ifit1", "Ifit2", "Ifit3",
    "Rsad2", "Ddx58", "Ifih1", "Stat1", "Stat2", "Irf9", "Ifi44", "Ifi27",
]
ISG_HUMAN = [
    "CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "IRF7",
    "CXCL10", "ISG15", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2",
    "IFIT3", "RSAD2", "DDX58", "IFIH1", "STAT1", "STAT2", "IRF9",
    "IFI44", "IFI27", "MX2",
]
RER_MOUSE = [
    "Rnaseh2a", "Rnaseh2b", "Rnaseh2c", "Rnaseh1", "Fen1", "Pold1",
    "Lig1", "Top1", "Tdp1", "Top2a",
]
RER_HUMAN = [
    "RNASEH2A", "RNASEH2B", "RNASEH2C", "RNASEH1", "FEN1", "POLD1",
    "LIG1", "TOP1", "TDP1", "TOP2A",
]
FTD_MOUSE = [
    "Mapt", "Grn", "C9orf72", "Tardbp", "Vcp", "Chmp2b", "Tbk1",
    "Spp1", "Tmem106b", "Gsk3b", "Cdk5", "Mark2",
]
FTD_HUMAN = [
    "MAPT", "GRN", "C9orf72", "TARDBP", "VCP", "CHMP2B", "TBK1",
    "SPP1", "TMEM106B", "GSK3B", "CDK5", "MARK2", "MARK3",
]
ONEC_HUMAN = [
    "MTHFR", "MTHFD1", "MTHFD2", "DHFR", "TYMS", "SHMT1", "SHMT2",
    "ATIC", "GART", "PPAT",
]
ENTREZ = {
    115004: "CGAS", 340061: "STING1", 29110: "TBK1", 3661: "IRF3",
    3665: "IRF7", 3627: "CXCL10", 9636: "ISG15", 4599: "MX1",
    4938: "OAS1", 4939: "OAS2", 4940: "OAS3", 3434: "IFIT1",
    3433: "IFIT2", 3437: "IFIT3", 91543: "RSAD2", 23586: "DDX58",
    64135: "IFIH1", 6772: "STAT1", 6773: "STAT2", 10379: "IRF9",
    10561: "IFI44", 3429: "IFI27", 4600: "MX2", 79621: "RNASEH2B",
    10535: "RNASEH2A", 84153: "RNASEH2C", 246243: "RNASEH1",
    2237: "FEN1", 5424: "POLD1", 3978: "LIG1", 7150: "TOP1",
    55775: "TDP1", 7153: "TOP2A", 4137: "MAPT", 2896: "GRN",
    203228: "C9orf72", 23435: "TARDBP", 7415: "VCP", 25978: "CHMP2B",
    6696: "SPP1", 54664: "TMEM106B", 2932: "GSK3B", 1020: "CDK5",
    2011: "MARK2", 4140: "MARK3", 4524: "MTHFR", 4522: "MTHFD1",
    10797: "MTHFD2", 1719: "DHFR", 7298: "TYMS", 6470: "SHMT1",
    6472: "SHMT2", 471: "ATIC", 2618: "GART", 5471: "PPAT",
}


def fnum(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def mean(xs):
    xs = [x for x in xs if x == x]
    return float(np.mean(xs)) if xs else float("nan")


def mwu(a, b):
    a = [x for x in a if x == x]
    b = [x for x in b if x == x]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    try:
        r = stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(r.statistic), float(r.pvalue)
    except Exception:
        return float("nan"), float("nan")


def log2fc(num, den):
    if num != num or den != den or num <= 0 or den <= 0:
        return float("nan")
    return math.log2(num / den)


def bedgraph_sum(path: Path) -> float:
    tot = 0.0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            p = line.split()
            if len(p) < 4:
                continue
            try:
                start, end, val = int(p[1]), int(p[2]), float(p[3])
            except ValueError:
                continue
            tot += val * max(end - start, 1)
    return tot


def test1_chrM():
    idx = list(csv.DictReader((PUB / "chrM/index.tsv").open(), delimiter="\t"))
    by = defaultdict(lambda: {"KOH": 0.0, "KCl": 0.0})
    for r in idx:
        p = PUB / "chrM" / r["file"]
        if not p.exists():
            continue
        key = (r["age_id"][0], r["age_id"], r["tissue"])
        by[key][r["chem"]] += bedgraph_sum(p)
    rows = []
    for (age, aid, tis), d in sorted(by.items()):
        koh, kcl = d["KOH"], d["KCl"]
        rows.append({
            "age": "old" if age == "o" else "young",
            "animal": aid,
            "tissue": tis,
            "koh": koh,
            "kcl": kcl,
            "koh_over_kcl": (koh / kcl) if kcl > 0 else float("nan"),
            "koh_minus_kcl": koh - kcl,
        })
    summary = {}
    for tis in ("brain", "heart", "kidney"):
        old = [r["koh_over_kcl"] for r in rows if r["tissue"] == tis and r["age"] == "old"]
        yng = [r["koh_over_kcl"] for r in rows if r["tissue"] == tis and r["age"] == "young"]
        u, p = mwu(old, yng)
        summary[tis] = {
            "old_mean_KOH_KCl": mean(old),
            "young_mean_KOH_KCl": mean(yng),
            "old_over_young": (mean(old) / mean(yng)) if mean(yng) else float("nan"),
            "mwu_U": u,
            "mwu_p": p,
            "n_old": len(old),
            "n_young": len(yng),
            "old": old,
            "young": yng,
        }
    return {"per_animal": rows, "by_tissue": summary}


def load_tsv_genes(path, gene_col, skip_hash=True):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        if skip_hash:
            peek = fh.readline()
            if not peek.startswith("#"):
                fh.seek(0)
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
        return rdr.fieldnames, rows


def group_means(rows, gene_col, groups):
    """groups: name -> list of column names. Return gene -> {group: mean}."""
    out = {}
    for r in rows:
        g = r[gene_col]
        out[g] = {name: mean([fnum(r[c]) for c in cols if c in r]) for name, cols in groups.items()}
    return out


def set_contrast(gene_means, genes, a, b):
    recs = []
    fcs = []
    for g in genes:
        if g not in gene_means:
            recs.append({"gene": g, "present": False})
            continue
        va, vb = gene_means[g].get(a), gene_means[g].get(b)
        fc = log2fc(va, vb) if va == va and vb == vb else float("nan")
        recs.append({"gene": g, "present": True, a: va, b: vb, "log2fc": fc})
        if fc == fc:
            fcs.append(fc)
    return {
        "n_present": sum(1 for r in recs if r.get("present")),
        "median_log2fc": float(np.median(fcs)) if fcs else float("nan"),
        "mean_log2fc": mean(fcs),
        "n_up_0.5": sum(1 for x in fcs if x >= 0.5),
        "n_down_0.5": sum(1 for x in fcs if x <= -0.5),
        "genes": recs,
    }


def test2_brain_h2():
    # cerebellum 9wk log2CPM already log2 — convert contrast as delta, not log2fc of CPM
    hdr, rows = load_tsv_genes(PUB / "GSE181159_9wk_log2CPM.txt.gz", "Gene")
    ctrl = [c for c in hdr if "CNTRL" in c]
    mut = [c for c in hdr if "Rnaseh2b" in c]
    out_sets = {}
    for label, genes in (("ISG", ISG_MOUSE), ("RER", RER_MOUSE), ("FTD_tau", FTD_MOUSE)):
        recs, deltas = [], []
        for g in genes:
            hits = [r for r in rows if r["Gene"] == g]
            if not hits:
                recs.append({"gene": g, "present": False})
                continue
            r = hits[0]
            c = mean([fnum(r[c]) for c in ctrl])
            m = mean([fnum(r[c]) for c in mut])
            d = m - c
            recs.append({"gene": g, "present": True, "WT": c, "Rnaseh2b": m, "delta_log2CPM": d})
            if d == d:
                deltas.append(d)
        out_sets[label] = {
            "n_present": sum(1 for r in recs if r.get("present")),
            "median_delta_log2CPM": float(np.median(deltas)) if deltas else float("nan"),
            "n_up_0.5": sum(1 for x in deltas if x >= 0.5),
            "n_down_0.5": sum(1 for x in deltas if x <= -0.5),
            "genes": recs,
        }

    # GSE193711 RNASEH2B-KO neurons + astrocytes (raw counts → CPM-ish via col sum)
    with gzip.open(PUB / "GSE193711_featureCounts.txt.gz", "rt") as fh:
        fh.readline()
        hdr = fh.readline().rstrip().split("\t")
        mat = list(csv.DictReader(fh, fieldnames=hdr, delimiter="\t"))
    colmap = {
        "neu_WT": [c for c in hdr if c.startswith("neurons_A6")],
        "neu_H2": [c for c in hdr if c.startswith("neurons_R19")],
        "neu_TX": [c for c in hdr if c.startswith("neurons_T13")],
        "ast_WT": [c for c in hdr if c.startswith("astro_A6")],
        "ast_H2": [c for c in hdr if c.startswith("astro_R19")],
        "ast_TX": [c for c in hdr if c.startswith("astro_T13")],
    }
    # library sizes
    lib = {k: 0.0 for k in colmap}
    for r in mat:
        for k, cols in colmap.items():
            lib[k] += sum(fnum(r[c]) for c in cols if c in r)
    gene_cpm = {}
    for r in mat:
        g = r["Geneid"]
        gene_cpm[g] = {}
        for k, cols in colmap.items():
            raw = mean([fnum(r[c]) for c in cols if c in r])
            gene_cpm[g][k] = (1e6 * raw / (lib[k] / max(len(cols), 1))) if lib[k] else float("nan")

    human = {}
    for tissue, a, b in (("neurons_H2_vs_WT", "neu_H2", "neu_WT"),
                         ("neurons_TREX1_vs_WT", "neu_TX", "neu_WT"),
                         ("astro_H2_vs_WT", "ast_H2", "ast_WT"),
                         ("astro_TREX1_vs_WT", "ast_TX", "ast_WT")):
        human[tissue] = {
            s: set_contrast(gene_cpm, genes, a, b)
            for s, genes in (("ISG", ISG_HUMAN), ("RER", RER_HUMAN), ("FTD_tau", FTD_HUMAN))
        }
    return {"GSE181159_cerebellum_9wk_Rnaseh2b": out_sets, "GSE193711_AGS_iPSC": human}


def test3_human_ribome():
    # 500 kb bins
    chr_sum = defaultdict(float)
    mapt_bin = 0.0
    n_bins = 0
    with gzip.open(PUB / "GSE297046_DNA_rez_500K_counts.tsv.gz", "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            n_bins += 1
            c = r["Chromosome"]
            chr_sum[c] += fnum(r["Count"])
            start, end = int(r["Start"]), int(r["End"])
            # hg38 MAPT ~ chr17:45.89-46.03 Mb
            if c.replace("chr", "") == "17" and start < 46_100_000 and end > 45_800_000:
                mapt_bin += fnum(r["Count"])
    tot = sum(chr_sum.values()) or 1.0
    mito_keys = [k for k in chr_sum if k.replace("chr", "").upper() in ("M", "MT")]
    mito = sum(chr_sum[k] for k in mito_keys)
    # HEK TPM: RNASEH2-related KO vs WT (column names)
    with gzip.open(PUB / "GSE297046_HEK_tpm.csv.gz", "rt") as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        ko = [c for c in fields if "KO" in c]
        wt = [c for c in fields if "WT" in c]
        tpm = {}
        for r in rdr:
            # Geneid is Ensembl; no symbol. Skip gene-set here unless we match later.
            tpm[r["Geneid"]] = {
                "KO": mean([fnum(r[c]) for c in ko]),
                "WT": mean([fnum(r[c]) for c in wt]),
            }
    return {
        "n_bins": n_bins,
        "total_counts": tot,
        "mito_fraction": mito / tot,
        "mito_chroms": mito_keys,
        "top_chrom_frac": sorted(
            ((k, v / tot) for k, v in chr_sum.items()), key=lambda x: -x[1]
        )[:8],
        "MAPT_chr17_45.8-46.1Mb_count": mapt_bin,
        "MAPT_bin_frac": mapt_bin / tot,
        "HEK_tpm_n_genes": len(tpm),
        "HEK_note": "TPM is RNA, not rNMP; KO vs WT columns present but gene IDs are Ensembl.",
        "HEK_cols_KO": ko,
        "HEK_cols_WT": wt,
    }


def drip_mapt():
    """Overlap DRIP peaks with hg38 MAPT (17:45894554-46028334) and TOP1 (20:41028816-41124476)."""
    loci = {
        "MAPT": ("17", 45_894_554, 46_028_334),
        "TOP1": ("20", 41_028_816, 41_124_476),
        "POLR2A": ("17", 7_483_983, 7_514_316),
        "RAD51": ("15", 40_694_753, 40_732_377),
    }
    files = {
        "siLuc": PUB / "GSE146970_siLuc.bed.gz",
        "siH2A": PUB / "GSE146970_siH2A.bed.gz",
        "H2A_ChIP": PUB / "GSE146970_H2A_chip.bed.gz",
    }

    def n_overlap(path, chrom, a, b):
        n = 0
        tot = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("track"):
                    continue
                p = line.split()
                if len(p) < 3:
                    continue
                tot += 1
                c = p[0].replace("chr", "")
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                if c == chrom and s < b and e > a:
                    n += 1
        return n, tot

    out = {}
    for lab, path in files.items():
        out[lab] = {}
        for name, (c, a, b) in loci.items():
            n, tot = n_overlap(path, c, a, b)
            out[lab][name] = {"peaks_overlapping": n, "total_peaks": tot}
    return out


def test5_drugs():
    # MTX: files Mo_{patient}_{A|B|C}.txt.gz  A=pre B=5d C=30d; gene_id=entrez
    mtx_dir = PUB / "gse279719"
    by_pt = defaultdict(dict)
    for p in sorted(mtx_dir.glob("*.gz")):
        m = re.search(r"Mo_(\d+)_([ABC])", p.name)
        if not m:
            continue
        pt, arm = m.group(1), m.group(2)
        raw = gzip.open(p, "rb").read()
        if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
            text = raw.decode("utf-16")
        else:
            text = raw.decode("utf-8", errors="replace")
        d = {}
        for line in text.splitlines()[1:]:
            a = line.split()
            if len(a) < 2:
                continue
            try:
                eid = int(float(a[0]))
            except ValueError:
                continue
            if eid in ENTREZ:
                d[ENTREZ[eid]] = fnum(a[1])
        by_pt[pt][arm] = d

    def paired(arm):
        sets = {}
        for label, genes in (("ISG", ISG_HUMAN), ("RER", RER_HUMAN),
                             ("FTD_tau", FTD_HUMAN), ("one_carbon", ONEC_HUMAN)):
            recs, fcs = [], []
            for g in genes:
                pre, post = [], []
                for pt, arms in by_pt.items():
                    if "A" in arms and arm in arms:
                        pre.append(arms["A"].get(g, float("nan")))
                        post.append(arms[arm].get(g, float("nan")))
                # paired log2fc per patient then median
                pfc = [log2fc(b, a) for a, b in zip(pre, post)]
                pfc = [x for x in pfc if x == x]
                recs.append({
                    "gene": g,
                    "n_pairs": len(pfc),
                    "median_log2fc": float(np.median(pfc)) if pfc else float("nan"),
                    "mean_pre": mean(pre),
                    "mean_post": mean(post),
                })
                fcs.extend(pfc)
            sets[label] = {
                "n_pairs_total": len(by_pt),
                "median_log2fc": float(np.median(fcs)) if fcs else float("nan"),
                "n_up_0.5": sum(1 for r in recs if (r["median_log2fc"] == r["median_log2fc"] and r["median_log2fc"] >= 0.5)),
                "n_down_0.5": sum(1 for r in recs if (r["median_log2fc"] == r["median_log2fc"] and r["median_log2fc"] <= -0.5)),
                "genes": recs,
            }
        return sets

    mtx = {"d5_vs_pre": paired("B"), "d30_vs_pre": paired("C"), "n_patients": len(by_pt)}

    # pemetrexed FPKM: L=LPS, PL=PMX+LPS, FL=FA+LPS, PFL=PMX+FA+LPS
    pem_dir = PUB / "gse159380"
    pem = defaultdict(dict)
    for p in pem_dir.glob("*.gz"):
        m = re.search(r"_(L|PL|FL|PFL)_([ABC])_fpkm", p.name)
        if not m:
            continue
        cond, rep = m.group(1), m.group(2)
        with gzip.open(p, "rt") as fh:
            next(fh)
            for line in fh:
                a = line.split()
                if len(a) < 2:
                    continue
                try:
                    eid = int(float(a[0]))
                except ValueError:
                    continue
                if eid in ENTREZ:
                    pem[(cond, rep)][ENTREZ[eid]] = fnum(a[1])

    def pem_contrast(a, b):
        sets = {}
        for label, genes in (("ISG", ISG_HUMAN), ("RER", RER_HUMAN),
                             ("FTD_tau", FTD_HUMAN), ("one_carbon", ONEC_HUMAN)):
            recs, fcs = [], []
            for g in genes:
                va = [pem[(a, r)].get(g, float("nan")) for r in "ABC"]
                vb = [pem[(b, r)].get(g, float("nan")) for r in "ABC"]
                fc = log2fc(mean(va), mean(vb))
                recs.append({"gene": g, a: mean(va), b: mean(vb), "log2fc": fc})
                if fc == fc:
                    fcs.append(fc)
            sets[label] = {
                "median_log2fc": float(np.median(fcs)) if fcs else float("nan"),
                "n_up_0.5": sum(1 for x in fcs if x >= 0.5),
                "n_down_0.5": sum(1 for x in fcs if x <= -0.5),
                "genes": recs,
            }
        return sets

    # GSE276623 RH2B OE: BC=RH2B ctrl, BD=RH2B dox, PC=parental ctrl, PD=parental dox
    with gzip.open(PUB / "GSE276623_Raw_counts.txt.gz", "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
    cols = {
        "RH2B_ctrl": ["BC1", "BC2", "BC3"],
        "RH2B_dox": ["BD1", "BD2", "BD3"],
        "par_ctrl": ["PC1", "PC2", "PC3"],
        "par_dox": ["PD1", "PD2", "PD3"],
    }
    lib = {k: 0.0 for k in cols}
    for r in rows:
        for k, cs in cols.items():
            lib[k] += sum(fnum(r[c]) for c in cs)
    gmean = {}
    for r in rows:
        g = r["SYMBOL"]
        gmean[g] = {
            k: 1e6 * mean([fnum(r[c]) for c in cs]) / (lib[k] / 3) if lib[k] else float("nan")
            for k, cs in cols.items()
        }
    rh2b = {
        s: set_contrast(gmean, genes, "RH2B_dox", "RH2B_ctrl")
        for s, genes in (("ISG", ISG_HUMAN), ("RER", RER_HUMAN), ("FTD_tau", FTD_HUMAN))
    }
    return {
        "GSE279719_MTX_monocytes": mtx,
        "GSE159380_PMX_macrophages": {
            "PMX_vs_LPS": pem_contrast("PL", "L"),
            "PMX_FA_vs_LPS": pem_contrast("PFL", "L"),
        },
        "GSE276623_RH2B_OE_HCT116": rh2b,
        "GSE276623_note": "Deposited RNA-seq is RNaseH2B overexpression ± dox, not HU treatment.",
    }


def verdicts(blob):
    t1 = blob["test1_chrM"]["by_tissue"]
    brain = t1["brain"]
    t2c = blob["test2_ISG"]["GSE181159_cerebellum_9wk_Rnaseh2b"]["ISG"]
    t2n = blob["test2_ISG"]["GSE193711_AGS_iPSC"]["neurons_H2_vs_WT"]["ISG"]
    t2a = blob["test2_ISG"]["GSE193711_AGS_iPSC"]["astro_H2_vs_WT"]["ISG"]
    t2nf = blob["test2_ISG"]["GSE193711_AGS_iPSC"]["neurons_H2_vs_WT"]["FTD_tau"]
    t3 = blob["test3_ribome"]
    t5m = blob["test5_drugs"]["GSE279719_MTX_monocytes"]["d30_vs_pre"]["ISG"]
    t5p = blob["test5_drugs"]["GSE159380_PMX_macrophages"]["PMX_vs_LPS"]["ISG"]
    t5r = blob["test5_drugs"]["GSE279719_MTX_monocytes"]["d30_vs_pre"]["RER"]
    return {
        "T1_aging_mtDNA_rNMP_brain": {
            "old_over_young_KOH_KCl": brain["old_over_young"],
            "p": brain["mwu_p"],
            "call": (
                "SUPPORTED_DIRECTION" if brain["old_over_young"] > 1.1 and (brain["mwu_p"] != brain["mwu_p"] or brain["mwu_p"] < 0.2)
                else "NOT_SUPPORTED_AT_N3"
                if brain["old_over_young"] <= 1.05
                else "WEAK_OR_UNDERPOWERED"
            ),
            "license": "GSE183589 WT mouse 7.5 vs 30 wk; n=3. Not human neurons. 30 wk is adult, not geriatric.",
        },
        "T2_H2_loss_ISG_brain": {
            "cerebellum_median_delta_log2CPM": t2c["median_delta_log2CPM"],
            "cerebellum_n_up": t2c["n_up_0.5"],
            "neuron_H2_median_log2fc": t2n["median_log2fc"],
            "astro_H2_median_log2fc": t2a["median_log2fc"],
            "call": (
                "SUPPORTED_IN_AGS_NEURAL"
                if (t2n["median_log2fc"] >= 0.3 or t2a["median_log2fc"] >= 0.3 or t2c["median_delta_log2CPM"] >= 0.3)
                else "NOT_SUPPORTED"
            ),
            "license": "RNase H2 loss in cerebellum / iPSC neurons+astrocytes. AGS, not FTD.",
        },
        "T3_nuclear_rNMP_MAPT": {
            "mito_fraction_500kb": t3["mito_fraction"],
            "MAPT_bin_frac": t3["MAPT_bin_frac"],
            "call": "BIN_TOO_COARSE_NO_GENE_BODY_MAP",
            "license": "500 kb bins cannot score long-gene bodies. MAPT occupancy is one bin.",
        },
        "T4_FTD_genes_after_H2_loss": {
            "neuron_H2_FTD_median_log2fc": t2nf["median_log2fc"],
            "neuron_H2_FTD_n_up": t2nf["n_up_0.5"],
            "call": (
                "EMPTY_OVERLAP"
                if (t2nf["median_log2fc"] != t2nf["median_log2fc"] or abs(t2nf["median_log2fc"]) < 0.25)
                else "SOME_MOVEMENT_NOT_A_PATH"
            ),
            "license": "MAPT/GRN/VCP etc. in RNASEH2B-KO neurons. Movement ≠ FTD mechanism.",
        },
        "T5_1C_drugs_ISG_RER": {
            "MTX_d30_ISG_median": t5m["median_log2fc"],
            "MTX_d30_RER_median": t5r["median_log2fc"],
            "PMX_ISG_median": t5p["median_log2fc"],
            "call": (
                "NO_ISG_INDUCTION"
                if (t5m["median_log2fc"] < 0.25 and t5p["median_log2fc"] < 0.25)
                else "MIXED"
            ),
            "license": "Patient monocytes on MTX; LPS macrophages + PMX. Not neurons. Not rNMP maps.",
        },
    }


def main():
    blob = {
        "test1_chrM": test1_chrM(),
        "test2_ISG": test2_brain_h2(),
        "test3_ribome": test3_human_ribome(),
        "test3b_DRIP_loci": drip_mapt(),
        "test5_drugs": test5_drugs(),
    }
    blob["verdicts"] = verdicts(blob)
    (OUT / "TEST_RESULTS.json").write_text(json.dumps(blob, indent=2, default=str))
    print(json.dumps(blob["verdicts"], indent=2, default=str))
    print("WROTE", OUT / "TEST_RESULTS.json")


if __name__ == "__main__":
    main()
