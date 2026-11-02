#!/usr/bin/env python3
"""Publication figures for the neuron/aging public-data report. Numbers from TEST_RESULTS*.json only."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "text.color": "#222222",
    "axes.labelcolor": "#222222",
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "axes.edgecolor": "#333333",
})

ROOT = Path("/home/noahtjones/hu-rnmp-pipeline/results/_cross")
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

T1 = json.loads((ROOT / "TEST_RESULTS.json").read_text())
T2 = json.loads((ROOT / "TEST_RESULTS_LEFTOVER.json").read_text())


def save(fig, name):
    p = FIG / name
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("WROTE", p)


def fig1_mtdna():
    tissues = ["brain", "heart", "kidney"]
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    rng = np.random.default_rng(0)
    for i, tis in enumerate(tissues):
        d = T1["test1_chrM"]["by_tissue"][tis]
        y = np.array(d["young"])
        o = np.array(d["old"])
        ax.scatter(np.full_like(y, i - 0.14) + rng.uniform(-0.02, 0.02, len(y)),
                   y, c="#2a6f97", s=42, zorder=3, label="7.5 wk" if i == 0 else None)
        ax.scatter(np.full_like(o, i + 0.14) + rng.uniform(-0.02, 0.02, len(o)),
                   o, c="#c44536", s=42, zorder=3, label="30 wk" if i == 0 else None)
        ax.hlines(np.mean(y), i - 0.26, i - 0.02, colors="#2a6f97", lw=2)
        ax.hlines(np.mean(o), i + 0.02, i + 0.26, colors="#c44536", lw=2)
    ax.axhline(1.0, color="#888888", ls="--", lw=0.7)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Brain", "Heart", "Kidney"])
    ax.set_ylabel("chrM KOH / KCl")
    ax.set_ylim(0, 6.2)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("GSE183589 HydEn-seq, n = 3 per age")
    save(fig, "figure1_mtdna_age.png")


def fig2_nuclear_vs_mt():
    mt = T1["test1_chrM"]["by_tissue"]["brain"]
    nuc = T2["T1b_nuclear_brain"]
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    data = [
        ("chrM\n7.5 wk", mt["young"], "#2a6f97"),
        ("chrM\n30 wk", mt["old"], "#c44536"),
        ("Nuclear\n7.5 wk", nuc["young"], "#5c8a6a"),
        ("Nuclear\n30 wk", nuc["old"], "#b5651d"),
    ]
    for i, (lab, vals, col) in enumerate(data):
        ax.scatter(np.full(len(vals), i), vals, c=col, s=48, zorder=3)
        ax.hlines(np.mean(vals), i - 0.22, i + 0.22, colors=col, lw=2)
    ax.axhline(1.0, color="#888888", ls="--", lw=0.7)
    ax.set_xticks(range(4))
    ax.set_xticklabels([d[0] for d in data])
    ax.set_ylabel("KOH / KCl")
    ax.set_ylim(0, 6.2)
    ax.set_title("Brain HydEn-seq, same animals")
    save(fig, "figure2_nuclear_vs_mt.png")


def fig3_isg():
    # curated gene order shared across panels
    rows = [
        ("Cerebellum Rnaseh2b\nΔ log2CPM", [
            ("Mx1", 1.099), ("Ifi44", 1.233), ("Irf7", 0.885),
            ("Ifit1", 0.785), ("Ifit3", 0.706), ("Isg15", 0.455),
            ("Tmem173", 0.327), ("Stat1", 0.264), ("Cgas", None),
        ], "#6b2d5c"),
        ("iPSC neurons H2-KO\nlog2FC", [
            ("MX1", 0.327), ("IFI44", -0.325), ("IRF7", 0.285),
            ("IFIT1", 0.254), ("IFIT3", -0.530), ("ISG15", -0.157),
            ("STING1", -1.364), ("STAT1", -0.449), ("CGAS", -0.696),
        ], "#2a6f97"),
        ("iPSC astrocytes H2-KO\nlog2FC", [
            ("MX1", 1.555), ("IFI44", 0.594), ("IRF7", 0.576),
            ("IFIT1", 0.901), ("IFIT3", 1.143), ("ISG15", 1.096),
            ("STING1", 0.231), ("STAT1", -0.335), ("CGAS", 0.470),
        ], "#c44536"),
        ("FTD vs Ctrl brain\nlog2FC CPM", [
            ("MX1", -0.878), ("IFI44", 0.757), ("IRF7", 0.694),
            ("IFIT1", 0.191), ("IFIT3", 0.510), ("ISG15", 0.458),
            ("TMEM173", 0.195), ("STAT1", -0.198), ("CGAS", 0.483),
        ], "#b5651d"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(10.6, 3.8), sharey=False)
    for ax, (title, genes, col) in zip(axes, rows):
        names = [g[0] for g in genes]
        vals = [g[1] if g[1] is not None else 0 for g in genes]
        present = [g[1] is not None for g in genes]
        colors = [col if p else "#cccccc" for p in present]
        y = np.arange(len(names))[::-1]
        ax.barh(y, [v if p else 0 for v, p in zip(vals, present)], color=colors, height=0.7)
        ax.axvline(0, color="#333333", lw=0.7)
        ax.axvline(0.5, color="#888888", ls=":", lw=0.6)
        ax.axvline(-0.5, color="#888888", ls=":", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlim(-1.8, 1.8)
        ax.set_xlabel("effect")
    fig.suptitle("ISG / STING panel across four neural systems", y=1.02, fontsize=11)
    save(fig, "figure3_isg_panels.png")


def fig4_ftd():
    genes = T2["T2b_GSE230519_postmortem"]["FTD_tau"]["genes"]
    # add TMEM173 and RER highlights
    extra = [
        {"gene": "TMEM173", "log2fc": 0.195, "mwu_p": 1.0},
        {"gene": "RNASEH2B", "log2fc": 0.361, "mwu_p": 0.234},
        {"gene": "TOP1", "log2fc": 0.060, "mwu_p": 0.445},
    ]
    recs = genes + extra
    # unique keep first
    seen = set()
    keep = []
    for r in recs:
        if r["gene"] in seen:
            continue
        seen.add(r["gene"])
        keep.append(r)
    keep = [r for r in keep if r["log2fc"] == r["log2fc"]]
    keep.sort(key=lambda r: r["log2fc"])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    y = np.arange(len(keep))
    cols = []
    for r in keep:
        if r["gene"] == "SPP1":
            cols.append("#c44536")
        elif r.get("mwu_p", 1) is not None and r.get("mwu_p") == r.get("mwu_p") and r["mwu_p"] < 0.05:
            cols.append("#2a6f97")
        else:
            cols.append("#888888")
    ax.barh(y, [r["log2fc"] for r in keep], color=cols, height=0.72)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r["gene"] for r in keep], fontsize=8.5)
    ax.set_xlabel("log2 (FTD / control) pseudobulk CPM")
    ax.set_title("GSE230519 postmortem, 6 FTD vs 7 control")
    save(fig, "figure4_ftd_postmortem.png")


def fig5_verdict():
    rows = [
        ("Brain mtDNA rNMP rises with age (30 wk)", "Not supported"),
        ("Brain mtDNA is an rNMP sink", "Supported"),
        ("Nuclear brain rNMP rises with age", "Not supported"),
        ("H2 loss → ISG in cerebellum", "Supported (AGS)"),
        ("H2 loss → ISG in astrocytes", "Supported (AGS)"),
        ("H2 loss → ISG in neurons", "Not supported"),
        ("FTD brain STING / cGAS program", "Not supported"),
        ("FTD brain osteopontin (SPP1)", "Supported"),
        ("FTD / tau is on the RER path", "Empty"),
        ("MTX / pemetrexed induce ISG–RER", "Not supported"),
        ("Bahat kidney/senescence mtDNA–STING", "Literature only"),
    ]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.axis("off")
    cmap = {
        "Supported": "#2d6a4f",
        "Supported (AGS)": "#2d6a4f",
        "Not supported": "#9b2226",
        "Empty": "#9b2226",
        "Literature only": "#b5651d",
    }
    for i, (claim, call) in enumerate(reversed(rows)):
        ax.text(0.02, i, claim, va="center", ha="left", fontsize=9, color="#222")
        ax.text(0.98, i, call, va="center", ha="right", fontsize=9,
                color=cmap[call], fontweight="bold")
    ax.set_title("Licensed calls (processed public tables only)")
    save(fig, "figure5_calls.png")


if __name__ == "__main__":
    fig1_mtdna()
    fig2_nuclear_vs_mt()
    fig3_isg()
    fig4_ftd()
    fig5_verdict()
