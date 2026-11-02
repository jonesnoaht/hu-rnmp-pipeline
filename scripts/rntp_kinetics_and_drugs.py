#!/usr/bin/env python3
"""Competitive rNTP vs dNTP insertion model + drug RNA-seq on RER-failure genes.

Kinetics from Nick McElhinny et al. PNAS 2010 (PMID 20194773) Table 1–2 and
reported stable-incorporation frequencies. Discrimination values for Pol ε
bases other than C are interpolated from the paper's ranking (ε lowest;
dC/rC = 500) and are tagged 'approx'. Drug RNA-seq licenses transcription
only, not rNMP maps.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

ROOT = Path("/home/noahtjones/hu-rnmp-pipeline")
PUB = ROOT / "data/public"
OUT = ROOT / "results/_cross"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "text.color": "#222222",
    "axes.labelcolor": "#222222",
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "axes.edgecolor": "#333333",
})

# Nick McElhinny 2010 Table 1 (µM, log-phase yeast)
DNTP = {"A": 16.0, "C": 14.0, "G": 12.0, "T": 30.0}
RNTP = {"A": 3000.0, "C": 500.0, "G": 700.0, "U": 1700.0}
PAIR = {"A": "A", "C": "C", "G": "G", "T": "U"}

# Discrimination D = (kcat/Km)_d / (kcat/Km)_r
# δ A/C from Table 2; δ G/T and α/ε from Fig. 1C text.
D_FACTOR = {
    "Pol δ": {"A": 1.6e4, "C": 1.3e4, "G": 1.0e4, "T": 1.7e6},
    "Pol ε": {"A": 2.0e3, "C": 5.0e2, "G": 2.0e3, "T": 5.0e3},  # C licensed; others approx
    "Pol α": {"A": 1.0e4, "C": 8.0e3, "G": 8.0e3, "T": 7.7e3},
}
D_SOURCE = {
    "Pol δ": {"A": "Table 2", "C": "Table 2", "G": "Fig 1C text", "T": "Fig 1C text"},
    "Pol ε": {"A": "approx", "C": "Fig 1C text", "G": "approx", "T": "approx"},
    "Pol α": {"A": "approx", "C": "approx", "G": "approx", "T": "Fig 1C vs δ"},
}

# Stable incorporation at physiological pools (paper Results)
STABLE_PER_DNMP = {"Pol ε": 1 / 1250, "Pol δ": 1 / 5000, "Pol α": 1 / 625}
# Yeast nuclear haploid genome, one S-phase
GENOME_BP = 1.2e7
# Rough polymerase share of newly made DNA (paper discussion)
POL_SHARE = {"Pol ε": 0.50, "Pol δ": 0.40, "Pol α": 0.10}

# Licensed 47-gene failed-excision set (NT≥10, HII≥10, ratio≥2)
TRUSTED = [
    "BOI1", "NAS2", "INP51", "CDC15", "FRE3", "POB3", "EMC1", "GSF2", "TNA1",
    "SGM1", "YJR098C", "SWC5", "GLK1", "IRA2", "MPH1", "VHS2", "ATR1", "UBP6",
    "SOK1", "HXT8", "ENT2", "CYC8", "SEN34", "FIP1", "AI1", "AI2", "tF(GAA)Q",
    "AI3", "YIL165C", "ATG13", "STP2", "TUP1", "BUB1", "BDP1", "ALG3",
    "tG(UCC)Q", "BUD14", "CTI6", "TAE1", "AI4", "DAL1", "GAT2", "YNL034W",
    "SCEI", "YCL068C", "DAL3", "YNR014W",
]
# Nuclear protein-coding only (drop mt ORFs, tRNAs, dubious) for drug tests
NUCLEAR_TRUSTED = [
    "BOI1", "NAS2", "INP51", "CDC15", "FRE3", "POB3", "EMC1", "GSF2", "TNA1",
    "SGM1", "SWC5", "GLK1", "IRA2", "MPH1", "VHS2", "ATR1", "UBP6", "SOK1",
    "HXT8", "ENT2", "CYC8", "SEN34", "FIP1", "ATG13", "STP2", "TUP1", "BUB1",
    "BDP1", "ALG3", "BUD14", "CTI6", "TAE1", "DAL1", "GAT2", "DAL3",
]
RER_YEAST = ["RNH201", "RNH202", "RNH203", "RNH1", "RAD27", "POL3", "CDC9", "TOP1", "TDP1"]

# High-confidence yeast→human (Alliance/SGD). Weak/no 1:1 omitted.
ORTHOLOGS = {
    "NAS2": [("PSMD9", "ENSG00000110801")],
    "INP51": [("SYNJ1", "ENSG00000159082"), ("SYNJ2", "ENSG00000078269")],
    "CDC15": [("STK4", "ENSG00000101109"), ("STK3", "ENSG00000104375")],
    "POB3": [("SSRP1", "ENSG00000149136")],
    "EMC1": [("EMC1", "ENSG00000127463")],
    "SWC5": [("CFDP1", "ENSG00000153774")],
    "GLK1": [("GCK", "ENSG00000106633"), ("HK1", "ENSG00000156515")],
    "IRA2": [("NF1", "ENSG00000196712")],
    "MPH1": [("FANCM", "ENSG00000187790")],
    "UBP6": [("USP14", "ENSG00000101557")],
    "ENT2": [("EPN1", "ENSG00000063245"), ("EPS15", "ENSG00000085832")],
    "CYC8": [("TLE1", "ENSG00000196781")],
    "TUP1": [("TLE1", "ENSG00000196781")],
    "SEN34": [("TSEN34", "ENSG00000170892")],
    "FIP1": [("FIP1L1", "ENSG00000145216")],
    "ATG13": [("ATG13", "ENSG00000175224")],
    "BUB1": [("BUB1", "ENSG00000169679")],
    "BDP1": [("BDP1", "ENSG00000145734")],
    "ALG3": [("ALG3", "ENSG00000214160")],
    "CTI6": [("PHF12", "ENSG00000109118")],
}
RER_HUMAN = {
    "RNASEH2A": "ENSG00000104889",
    "RNASEH2B": "ENSG00000136147",
    "RNASEH2C": "ENSG00000172922",
    "RNASEH1": "ENSG00000167525",
    "FEN1": "ENSG00000168496",
    "POLD1": "ENSG00000062822",
    "LIG1": "ENSG00000105486",
    "TOP1": "ENSG00000198900",
    "TDP1": "ENSG00000042088",
}
ISG_HUMAN = {
    "CGAS": "ENSG00000164430",
    "STING1": "ENSG00000184584",
    "TBK1": "ENSG00000183735",
    "IRF3": "ENSG00000126456",
    "IRF7": "ENSG00000185507",
    "CXCL10": "ENSG00000169245",
    "ISG15": "ENSG00000187608",
    "MX1": "ENSG00000157601",
    "OAS1": "ENSG00000089127",
    "OAS2": "ENSG00000111335",
    "OAS3": "ENSG00000111331",
    "IFIT1": "ENSG00000185745",
    "IFIT2": "ENSG00000119922",
    "IFIT3": "ENSG00000119917",
    "RSAD2": "ENSG00000134321",
    "DDX58": "ENSG00000107201",
    "IFIH1": "ENSG00000115267",
    "STAT1": "ENSG00000115415",
    "STAT2": "ENSG00000170581",
    "IRF9": "ENSG00000213928",
    "IFI44": "ENSG00000137965",
    "IFI27": "ENSG00000165949",
    "MX2": "ENSG00000183486",
}


def p_insert(D, r, d):
    if d <= 0:
        return 1.0
    return 1.0 / (1.0 + D * (d / r))


def insertions_per_sphase(d_scale=1.0, shares=None):
    """Competitive first-insertion load using D and scaled dNTP pools."""
    shares = shares or POL_SHARE
    total = 0.0
    by_pol = {}
    for pol, sh in shares.items():
        pr = []
        for base in ("A", "C", "G", "T"):
            pr.append(p_insert(D_FACTOR[pol][base], RNTP[PAIR[base]], DNTP[base] * d_scale))
        mean_p = float(np.mean(pr))
        n = GENOME_BP * sh * mean_p
        by_pol[pol] = {"mean_P": mean_p, "insertions": n}
        total += n
    return total, by_pol


def remaining(insertions, repair_frac):
    return insertions * (1.0 - repair_frac)


def fig_kinetics():
    # K1: P vs ratio
    ratios = np.logspace(0, 4, 80)
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    cols = {"Pol ε": "#c44536", "Pol δ": "#2a6f97", "Pol α": "#5c8a6a"}
    for pol, Ds in D_FACTOR.items():
        # show C (best-licensed ε number) and mean of four
        pC = [p_insert(Ds["C"], r, 1.0) for r in ratios]
        ax.plot(ratios, pC, color=cols[pol], lw=2, label=f"{pol} (C)")
    for base, r0, d0 in [("A", 3000, 16), ("C", 500, 14), ("G", 700, 12), ("T", 1700, 30)]:
        ax.axvline(r0 / d0, color="#bbbbbb", lw=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("[rNTP] / [dNTP]")
    ax.set_ylabel("P(insert rNMP | complementary)")
    ax.set_title("Competitive insertion")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "figureK1_competition.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # K2: HU-like dNTP depletion
    folds = np.array([1, 1.5, 2, 3, 5, 8, 10, 15, 20], dtype=float)
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    paper = 0.0
    for pol, sh in POL_SHARE.items():
        paper += GENOME_BP * sh * STABLE_PER_DNMP[pol]
    ax.plot(folds, paper * folds, color="#c44536", marker="o", lw=2, label="licensed stable × fold")
    ax.set_xlabel("Fold dNTP depletion (linear HU / RNR-block model)")
    ax.set_ylabel("rNMP insertions per S phase")
    ax.set_title("Scale the PNAS 2010 stable frequencies, not raw D")
    ax.legend(frameon=False)
    fig.savefig(FIG / "figureK2_hu_depletion.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # K3: remaining vs repair at WT and 10× HU using licensed stable load
    repair = np.linspace(0.0, 0.999, 80)
    wt, hu = paper, paper * 10
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(1 - repair, remaining(wt, repair), color="#2a6f97", lw=2, label="WT pools (~7,680)")
    ax.plot(1 - repair, remaining(hu, repair), color="#c44536", lw=2, label="10× dNTP drop")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Unrepaired fraction (1 − RER efficiency)")
    ax.set_ylabel("rNMPs remaining after S phase")
    ax.set_title("Pool stress × failed excision")
    ax.legend(frameon=False)
    fig.savefig(FIG / "figureK3_repair.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return


def model_table():
    rows = []
    for fold in (1, 2, 5, 10, 20):
        tot, by = insertions_per_sphase(1.0 / fold)
        rec = {"dNTP_fold_drop": fold, "insertions": tot}
        for pol, d in by.items():
            rec[f"{pol}_P"] = d["mean_P"]
            rec[f"{pol}_n"] = d["insertions"]
        rec["remain_RER_null"] = remaining(tot, 0.0)
        rec["remain_RER_99"] = remaining(tot, 0.99)
        rec["remain_RER_99.9"] = remaining(tot, 0.999)
        rows.append(rec)
    # also paper-stable frequencies (not D-derived)
    paper = 0.0
    for pol, sh in POL_SHARE.items():
        paper += GENOME_BP * sh * STABLE_PER_DNMP[pol]
    return {
        "pools_uM": {"dNTP": DNTP, "rNTP": RNTP},
        "ratios": {b: RNTP[PAIR[b]] / DNTP[b] for b in DNTP},
        "D": D_FACTOR,
        "D_source": D_SOURCE,
        "paper_stable_insertions_per_S": paper,
        "competitive_grid": rows,
        "license": (
            "D and pools: PMID 20194773. Competitive P = 1/(1+D[d]/[r]). "
            "Pol ε A/G/T D are approximate. HU is a uniform dNTP scale, not a measured pool table. "
            "MTX/pemetrexed are uracil-in-DNA, not this rNTP competition."
        ),
    }


def lfc(a, b):
    if a is None or b is None or a != a or b != b or a <= 0 or b <= 0:
        return float("nan")
    return math.log2(a / b)


def mwu(a, b):
    a = [x for x in a if x == x]
    b = [x for x in b if x == x]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def score_set(gene_fc: dict, genes):
    vals = [gene_fc[g] for g in genes if g in gene_fc and gene_fc[g] == gene_fc[g]]
    if not vals:
        return {"n": 0, "median": float("nan"), "n_up_0.5": 0, "n_down_0.5": 0, "genes": []}
    return {
        "n": len(vals),
        "median": float(np.median(vals)),
        "mean": float(np.mean(vals)),
        "n_up_0.5": int(sum(1 for v in vals if v >= 0.5)),
        "n_down_0.5": int(sum(1 for v in vals if v <= -0.5)),
        "genes": sorted(
            [{"gene": g, "log2fc": gene_fc[g]} for g in genes if g in gene_fc and gene_fc[g] == gene_fc[g]],
            key=lambda r: -abs(r["log2fc"]),
        ),
    }


def yeast_gse83454():
    """WT FPKM, t0 vs 1/2/4 h HU. Two reps."""
    d = PUB / "gse83454"
    groups = {
        "t0": ["GSM2203512_ctrl_t0_rep1.txt.gz", "GSM2203513_ctrl_t0_rep2.txt.gz"],
        "h1": ["GSM2203516_ctrl_1hHU_rep1.txt.gz", "GSM2203517_ctrl_1hHU_rep2.txt.gz"],
        "h2": ["GSM2203520_ctrl_2hHU_rep1.txt.gz", "GSM2203521_ctrl_2hHU_rep2.txt.gz"],
        "h4": ["GSM2203524_ctrl_4hHU_rep1.txt.gz", "GSM2203525_ctrl_4hHU_rep2.txt.gz"],
    }
    tables = {}
    for g, files in groups.items():
        acc = defaultdict(list)
        for fn in files:
            with gzip.open(d / fn, "rt") as fh:
                r = csv.DictReader(fh, delimiter="\t")
                for row in r:
                    gid = row["gene_id"]
                    try:
                        acc[gid].append(float(row["FPKM"]))
                    except (ValueError, KeyError):
                        continue
        tables[g] = {k: float(np.mean(v)) for k, v in acc.items()}

    def fc_map(treat):
        out = {}
        for g, tv in tables[treat].items():
            out[g] = lfc(tv, tables["t0"].get(g, float("nan")))
        return out

    blob = {}
    for treat in ("h1", "h2", "h4"):
        fc = fc_map(treat)
        # also by symbol via short name — cufflinks gene_id is systematic
        blob[treat] = {
            "trusted_sys": score_set(fc, NUCLEAR_TRUSTED),
            "rer": score_set(fc, RER_YEAST),
            "genome_median": float(np.nanmedian(list(fc.values()))),
            "n_genome": int(sum(1 for v in fc.values() if v == v)),
        }
        # map systematic via TRUSTED symbols: FPKM table uses Y-names; also gene_short_name
        # rebuild symbol fc
        sym = {}
        for fn in groups[treat]:
            break
        # read one treat + t0 with short names
        # simpler: parse t0 and treat again for short_name
        # already have systematic. Build symbol from a scan of t0 file
        pass

    # symbol map from first file
    sys_to_sym = {}
    with gzip.open(d / groups["t0"][0], "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            sn = row.get("gene_short_name") or ""
            if sn and sn != "-":
                sys_to_sym[row["gene_id"]] = sn
            else:
                sys_to_sym[row["gene_id"]] = row["gene_id"]

    for treat in ("h1", "h2", "h4"):
        fc_sys = fc_map(treat)
        fc_sym = {sys_to_sym.get(k, k): v for k, v in fc_sys.items()}
        # also keep systematic keys
        fc_both = dict(fc_sys)
        fc_both.update(fc_sym)
        blob[treat]["trusted"] = score_set(fc_both, NUCLEAR_TRUSTED)
        blob[treat]["rer"] = score_set(fc_both, RER_YEAST)
        blob[treat]["top_trusted"] = score_set(fc_both, NUCLEAR_TRUSTED)["genes"][:15]
    return blob


def yeast_gse63516():
    """log2 RPKM means. Contrast WT-HU vs WT-G1 (cell-cycle confounded) and vs WT-R60."""
    p = PUB / "gse63516/GSE63516_RNAseq_allmRNA_feature_report.txt.gz"
    rows = []
    with gzip.open(p, "rt", errors="replace") as fh:
        next(fh)  # title
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            feat = row.get("Feature") or ""
            try:
                rec = {
                    "id": feat,
                    "WT_G1": float(row["Mean WT-G1"]),
                    "WT_HU": float(row["Mean WT-HU"]),
                    "WT_R60": float(row["Mean WT-R60"]),
                }
            except (KeyError, ValueError):
                continue
            rows.append(rec)
    # table is already log2 RPKM — delta is log2FC if values are log2
    # "Log2 transformed" in header. delta HU-G1 is log2FC.
    fc_g1, fc_r60 = {}, {}
    for rec in rows:
        fc_g1[rec["id"]] = rec["WT_HU"] - rec["WT_G1"]
        fc_r60[rec["id"]] = rec["WT_HU"] - rec["WT_R60"]
    return {
        "note": "Values are log2 RPKM. HU vs G1 is S-arrest vs G1, not a clean drug contrast. HU vs R60 is on-drug vs 60 min recovery.",
        "HU_vs_G1_trusted": score_set(fc_g1, NUCLEAR_TRUSTED),
        "HU_vs_R60_trusted": score_set(fc_r60, NUCLEAR_TRUSTED),
        "HU_vs_G1_rer": score_set(fc_g1, RER_YEAST),
        "HU_vs_R60_rer": score_set(fc_r60, RER_YEAST),
        "HU_vs_G1_genome_median": float(np.nanmedian(list(fc_g1.values()))),
        "HU_vs_R60_genome_median": float(np.nanmedian(list(fc_r60.values()))),
    }


def human_gse229868():
    """HCT116 total RNA HU vs NT, 3 reps. Ensembl IDs."""
    p = PUB / "gse229868/GSE229868_HUtreatment.rnaseq.counts_12samples.txt.gz"
    want = {}
    for sym, ens in list(RER_HUMAN.items()) + list(ISG_HUMAN.items()):
        want[ens] = ("RER" if sym in RER_HUMAN else "ISG", sym)
    for y, pairs in ORTHOLOGS.items():
        for sym, ens in pairs:
            want[ens] = ("ORTH", f"{y}:{sym}")
    hu = ["Total_HU1", "Total_HU2", "Total_HU3"]
    nt = ["Total_NT1", "Total_NT2", "Total_NT3"]
    recs = {}
    lib = defaultdict(float)
    with gzip.open(p, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # first col blank
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            a = line.rstrip("\n").split("\t")
            gid = a[0].split(".")[0]
            try:
                vals = {h: float(a[idx[h]]) for h in hu + nt}
            except (ValueError, KeyError, IndexError):
                continue
            for h, v in vals.items():
                lib[h] += v
            if gid in want:
                recs[gid] = vals
    fc = {}
    grouped = {"ORTH": {}, "RER": {}, "ISG": {}}
    for ens, vals in recs.items():
        kind, name = want[ens]
        hu_cpm = [1e6 * vals[h] / lib[h] if lib[h] else float("nan") for h in hu]
        nt_cpm = [1e6 * vals[h] / lib[h] if lib[h] else float("nan") for h in nt]
        grouped[kind][name] = lfc(float(np.nanmean(hu_cpm)), float(np.nanmean(nt_cpm)))
        fc[name] = grouped[kind][name]
    return {
        "ORTH": score_set(grouped["ORTH"], list(grouped["ORTH"])),
        "RER": score_set(grouped["RER"], list(grouped["RER"])),
        "ISG": score_set(grouped["ISG"], list(grouped["ISG"])),
        "n_wanted": len(want),
        "n_found": len(recs),
    }


def human_gse274294():
    """A2058 siCTRL total RNA HU vs NT. Description order maps D029T*."""
    # titles in series-matrix order ↔ D029T11..T34
    titles = [
        "totalRNA_siCTRL_NT_REP1", "totalRNA_siCTRL_HU_REP1",
        "totalRNA_siLSM14A_NT_REP1", "totalRNA_siLSM14A_HU_REP1",
        "totalRNA_siCTRL_NT_REP2", "totalRNA_siCTRL_HU_REP2",
        "totalRNA_siLSM14A_NT_REP2", "totalRNA_siLSM14A_HU_REP2",
        "totalRNA_siCTRL_NT_REP3", "totalRNA_siCTRL_HU_REP3",
        "totalRNA_siLSM14A_NT_REP3", "totalRNA_siLSM14A_HU_REP3",
        "polysomalRNA_siCTRL_NT_REP1", "polysomalRNA_siCTRL_HU_REP1",
        "polysomalRNA_siLSM14A_NT_REP1", "polysomalRNA_siLSM14A_HU_REP1",
        "polysomalRNA_siCTRL_NT_REP2", "polysomalRNA_siCTRL_HU_REP2",
        "polysomalRNA_siLSM14A_NT_REP2", "polysomalRNA_siLSM14A_HU_REP2",
        "polysomalRNA_siCTRL_NT_REP3", "polysomalRNA_siCTRL_HU_REP3",
        "polysomalRNA_siLSM14A_NT_REP3", "polysomalRNA_siLSM14A_HU_REP3",
    ]
    desc = [f"D029T{i}" for i in range(11, 35)]
    title_of = dict(zip(desc, titles))
    nt_cols = [d for d, t in title_of.items() if t.startswith("totalRNA_siCTRL_NT")]
    hu_cols = [d for d, t in title_of.items() if t.startswith("totalRNA_siCTRL_HU")]
    want = {}
    for sym, ens in list(RER_HUMAN.items()) + list(ISG_HUMAN.items()):
        want[ens] = ("RER" if sym in RER_HUMAN else "ISG", sym)
    for y, pairs in ORTHOLOGS.items():
        for sym, ens in pairs:
            want[ens] = ("ORTH", f"{y}:{sym}")
    recs = {}
    lib = defaultdict(float)
    p = PUB / "gse274294/GSE274294_tablecounts_raw.csv.gz"
    with gzip.open(p, "rt") as fh:
        r = csv.DictReader(fh)
        for row in r:
            gid = row["geneID"].split(".")[0]
            for c in nt_cols + hu_cols:
                try:
                    v = float(row[c])
                except (KeyError, ValueError):
                    continue
                lib[c] += v
                if gid in want:
                    recs.setdefault(gid, {})[c] = v
    grouped = {"ORTH": {}, "RER": {}, "ISG": {}}
    for ens, vals in recs.items():
        kind, name = want[ens]
        hu_cpm = [1e6 * vals[c] / lib[c] for c in hu_cols if c in vals and lib[c]]
        nt_cpm = [1e6 * vals[c] / lib[c] for c in nt_cols if c in vals and lib[c]]
        if hu_cpm and nt_cpm:
            grouped[kind][name] = lfc(float(np.mean(hu_cpm)), float(np.mean(nt_cpm)))
    return {
        "nt_cols": nt_cols,
        "hu_cols": hu_cols,
        "ORTH": score_set(grouped["ORTH"], list(grouped["ORTH"])),
        "RER": score_set(grouped["RER"], list(grouped["RER"])),
        "ISG": score_set(grouped["ISG"], list(grouped["ISG"])),
    }


def fig_drugs(yeast, human_hu, human_mel, mtx_note):
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.8))
    # yeast 4h HU trusted genes
    ygenes = yeast["h4"]["trusted"]["genes"]
    names = [g["gene"] for g in ygenes]
    vals = [g["log2fc"] for g in ygenes]
    y = np.arange(len(names))[::-1]
    cols = ["#c44536" if v >= 0.5 else "#2a6f97" if v <= -0.5 else "#888888" for v in vals]
    axes[0].barh(y, vals, color=cols, height=0.7)
    axes[0].axvline(0, color="#333", lw=0.7)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(names, fontsize=6.5)
    axes[0].set_title("Yeast WT 4 h HU\nGSE83454 nuclear 47-set")
    axes[0].set_xlabel("log2FC vs t0")

    def panel(ax, blob, title):
        chunks = []
        for lab, key in [("ortholog", "ORTH"), ("RER", "RER"), ("ISG", "ISG")]:
            if key in blob and blob[key]["n"]:
                chunks.append((lab, blob[key]["median"], blob[key]["n"]))
        if not chunks:
            ax.text(0.5, 0.5, "no genes", ha="center")
            return
        xs = np.arange(len(chunks))
        ax.bar(xs, [c[1] for c in chunks], color=["#b5651d", "#2a6f97", "#6b2d5c"])
        ax.axhline(0, color="#333", lw=0.7)
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{c[0]}\nn={c[2]}" for c in chunks], fontsize=8)
        ax.set_ylabel("median log2FC")
        ax.set_title(title)

    panel(axes[1], human_hu, "HCT116 HU\nGSE229868 total RNA")
    panel(axes[2], human_mel, "A2058 HU\nGSE274294 siCTRL total")
    fig.tight_layout()
    fig.savefig(FIG / "figureK4_drug_sets.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    fig_kinetics()
    model = model_table()
    yeast = yeast_gse83454()
    y635 = yeast_gse63516()
    h229 = human_gse229868()
    h274 = human_gse274294()
    fig_drugs(yeast, h229, h274, None)
    blob = {
        "model": model,
        "yeast_GSE83454": yeast,
        "yeast_GSE63516": y635,
        "human_GSE229868_HCT116_HU": h229,
        "human_GSE274294_A2058_HU": h274,
        "orthologs_used": ORTHOLOGS,
        "verdicts": {
            "kinetics": (
                "At published yeast pools, competitive D predicts rNMP insertion "
                "that rises steeply when dNTP falls. RER-null leaves that load in place. "
                "MTX is not this competition."
            ),
            "yeast_HU_47": {
                "h4_median": yeast["h4"]["trusted"]["median"],
                "h4_n": yeast["h4"]["trusted"]["n"],
                "h4_up": yeast["h4"]["trusted"]["n_up_0.5"],
                "h4_down": yeast["h4"]["trusted"]["n_down_0.5"],
                "call": (
                    "NOT_A_COHERENT_PROGRAM"
                    if abs(yeast["h4"]["trusted"]["median"]) < 0.3
                    else "SHIFTED"
                ),
            },
            "human_HU_orthologs": {
                "HCT116_median": h229["ORTH"]["median"],
                "A2058_median": h274["ORTH"]["median"],
                "HCT116_ISG": h229["ISG"]["median"],
                "A2058_ISG": h274["ISG"]["median"],
            },
        },
    }
    (OUT / "KINETICS_AND_DRUGS.json").write_text(json.dumps(blob, indent=2, default=str))
    print(json.dumps(blob["verdicts"], indent=2, default=str))
    print("yeast 4h trusted", yeast["h4"]["trusted"])
    print("yeast 4h rer", yeast["h4"]["rer"])
    print("63516 HU vs R60", y635["HU_vs_R60_trusted"])
    print("HCT116", h229)
    print("A2058", {k: h274[k] for k in ("ORTH", "RER", "ISG")})
    print("WT insertions", model["competitive_grid"][0]["insertions"],
          "paper-stable", model["paper_stable_insertions_per_S"])


if __name__ == "__main__":
    main()
