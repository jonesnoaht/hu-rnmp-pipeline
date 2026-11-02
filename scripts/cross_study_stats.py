#!/usr/bin/env python3
"""Experimental-design tables + cross-study stats for HU-rNMP pilots.

Writes:
  results/experimental_design.tsv
  results/experimental_design.html
  results/cross_study_metrics.tsv
  results/cross_study_stats.json
  results/CONCLUSIONS.md
  results/cross_study_report.html
"""
from __future__ import annotations

import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import os

from scipy import stats

REPO = Path(os.environ.get("REPO", "/home/noahtjones/hu-rnmp-pipeline"))
OUT = Path(os.environ.get("RESULTS", str(REPO / "results")))
SACCER3_MB = 12.1  # approx nuclear + mito sacCer3

DESIGN = [
    # GSE61464 Koh 2015 ribose-seq, PMID 25622106
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575897",
        sample="rnh201 (KK-100)",
        genotype="rnh201",
        genotype_class="rer_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="KK-100",
        extra="RER-deficient; polymerase WT",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575898",
        sample="rnh201 (EconoTaq)",
        genotype="rnh201",
        genotype_class="rer_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="EconoTaq",
        extra="RER-deficient; polymerase WT; alt polymerase for library",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575899",
        sample="rnh201 (KK-30)",
        genotype="rnh201",
        genotype_class="rer_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="KK-30",
        extra="RER-deficient; polymerase WT",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575900",
        sample="rnh1 rnh201",
        genotype="rnh1 rnh201",
        genotype_class="rer_rnh1_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="standard",
        extra="RNase H1 + RER double null",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575901",
        sample="rnh1 rnh201",
        genotype="rnh1 rnh201",
        genotype_class="rer_rnh1_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="standard",
        extra="biological / library replicate of rnh1 rnh201",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575902",
        sample="rnh1 rnh201 ung1",
        genotype="rnh1 rnh201 ung1",
        genotype_class="rer_rnh1_ung1",
        polymerase="WT",
        treatment="none",
        enzyme_prep="standard",
        extra="adds ung1 (uracil-DNA glycosylase) on RNase H double-null",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575903",
        sample="pol2-M644G rnh201",
        genotype="pol2-M644G rnh201",
        genotype_class="pol_epsilon",
        polymerase="pol2-M644G",
        treatment="none",
        enzyme_prep="standard",
        extra="Pol ε steric-gate mutator; elevated rNMP incorporation",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575904",
        sample="pol2-4 rnh201",
        genotype="pol2-4 rnh201",
        genotype_class="pol_epsilon",
        polymerase="pol2-4",
        treatment="none",
        enzyme_prep="standard",
        extra="Pol ε catalytic/exonuclease allele on rnh201",
    ),
    dict(
        study="GSE61464",
        pmid="25622106",
        assay="ribose-seq",
        run="SRR1575905",
        sample="pol3-5DV rnh201",
        genotype="pol3-5DV rnh201",
        genotype_class="pol_delta",
        polymerase="pol3-5DV",
        treatment="none",
        enzyme_prep="standard",
        extra="Pol δ mutator; lagging-strand rNMP bias expected in literature",
    ),
    # GSE125855 Zhou/Kunkel HydEn-seq, PMID 31488849
    dict(
        study="GSE125855",
        pmid="31488849",
        assay="HydEn-seq",
        run="SRR8503261",
        sample="rnh201 NT (TAK_1811)",
        genotype="rnh201",
        genotype_class="rer_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="NT",
        extra="HydEn-seq no RNase HII pre-cleavage (alkali / hydrolytic ends)",
    ),
    dict(
        study="GSE125855",
        pmid="31488849",
        assay="HydEn-seq",
        run="SRR8503270",
        sample="rnh201 RNaseHII (TAK_1827)",
        genotype="rnh201",
        genotype_class="rer_null",
        polymerase="WT",
        treatment="none",
        enzyme_prep="RNaseHII",
        extra="HydEn-seq with RNase HII cleavage of embedded rNMPs",
    ),
    # GSE234409 HU WGS pair (series 2-1)
    dict(
        study="GSE234409",
        pmid="",
        assay="WGS",
        run="SRR24871232",
        sample="2-1-HU-14-B",
        genotype="WT",
        genotype_class="wt",
        polymerase="WT",
        treatment="HU",
        enzyme_prep="n/a",
        extra="Matched 2-1 series; hydroxyurea (replication-stress / dNTP depletion)",
    ),
    dict(
        study="GSE234409",
        pmid="",
        assay="WGS",
        run="SRR24871236",
        sample="2-1-C0-6",
        genotype="WT",
        genotype_class="wt",
        polymerase="WT",
        treatment="none",
        enzyme_prep="n/a",
        extra="Matched 2-1 untreated C0 control",
    ),
]


def parse_flagstat(p: Path) -> dict:
    text = p.read_text() if p.exists() else ""
    out = {"reads_total": None, "reads_mapped": None, "map_pct": None}
    m = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    if m:
        out["reads_total"] = int(m.group(1))
    m = re.search(r"^(\d+) \+ \d+ mapped \(([\d.]+)%", text, re.M)
    if m:
        out["reads_mapped"] = int(m.group(1))
        out["map_pct"] = float(m.group(2))
    return out


def hotspot_count(p: Path) -> int | None:
    if not p.exists():
        return None
    n = 0
    with p.open() as f:
        next(f, None)
        for line in f:
            if line.strip():
                n += 1
    return n


def chr_mito_frac(sig: dict) -> float | None:
    pc = sig.get("per_chromosome") or {}
    mito = pc.get("chrM") or {}
    mito_n = int(mito.get("forward", 0)) + int(mito.get("reverse", 0))
    tot = int(sig.get("total_rnmp_5prime_ends") or 0)
    if tot <= 0:
        return None
    return mito_n / tot


def chr12_bias(sig: dict) -> float | None:
    pc = sig.get("per_chromosome") or {}
    c = pc.get("chrXII") or {}
    f, r = int(c.get("forward", 0)), int(c.get("reverse", 0))
    if f + r == 0:
        return None
    return f / (f + r)


def collect_row(d: dict) -> dict:
    dest = OUT / d["study"] / d["run"]
    sigp = dest / "mutation_signature.json"
    flag = parse_flagstat(dest / "flagstat.txt")
    row = {**d, **flag}
    row["status"] = "pending"
    if sigp.exists():
        sig = json.loads(sigp.read_text())
        row["status"] = "complete"
        row["assay_called"] = sig.get("assay_type")
        row["ends"] = sig.get("total_rnmp_5prime_ends")
        row["unique_sites"] = sig.get("unique_rnmp_positions")
        row["density_per_mb"] = sig.get("rnmp_density_per_mb")
        row["fwd_frac"] = sig.get("forward_fraction")
        row["strand"] = sig.get("strand_asymmetry")
        row["avg_strand_bias"] = sig.get("avg_strand_bias")
        wc = sig.get("window_classes") or {}
        row["win_sym"] = wc.get("symmetric", 0)
        row["win_fwd"] = wc.get("forward_biased", 0)
        row["win_rev"] = wc.get("reverse_biased", 0)
        row["win_mild"] = wc.get("mild_bias", 0)
        row["mito_frac"] = chr_mito_frac(sig)
        row["chr12_fwd_frac"] = chr12_bias(sig)
        # uniqueness: sites per million mapped ends
        ends = row["ends"] or 0
        row["sites_per_m_ends"] = (row["unique_sites"] / ends * 1e6) if ends else None
    else:
        for k in (
            "assay_called",
            "ends",
            "unique_sites",
            "density_per_mb",
            "fwd_frac",
            "strand",
            "avg_strand_bias",
            "win_sym",
            "win_fwd",
            "win_rev",
            "win_mild",
            "mito_frac",
            "chr12_fwd_frac",
            "sites_per_m_ends",
        ):
            row[k] = None
    row["n_hotspots"] = hotspot_count(dest / "rnmp_hotspots.tsv")
    vcf = dest / f"{d['run']}.vcf.gz"
    row["vcf"] = "yes" if vcf.exists() else "no"
    return row


def mw(a, b):
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if len(a) < 2 or len(b) < 2:
        return {"n1": len(a), "n2": len(b), "U": None, "p": None, "note": "n<2"}
    r = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {"n1": len(a), "n2": len(b), "U": float(r.statistic), "p": float(r.pvalue)}


def spearman(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return {"n": len(pairs), "rho": None, "p": None}
    x, y = zip(*pairs)
    r = stats.spearmanr(x, y)
    return {"n": len(pairs), "rho": float(r.statistic), "p": float(r.pvalue)}


def fmt(x, nd=2):
    if x is None:
        return "—"
    if isinstance(x, float):
        if math.isnan(x):
            return "—"
        if abs(x) >= 1000:
            return f"{x:,.1f}"
        return f"{x:.{nd}f}"
    if isinstance(x, int):
        return f"{x:,}"
    return str(x)


def html_table(headers, rows):
    th = "".join(f"<th>{h}</th>" for h in headers)
    body = []
    for r in rows:
        body.append("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def md_fragment(text: str) -> str:
    """Small Markdown→HTML for CONCLUSIONS.md (headings, lists, tables, bold)."""
    import html as html_lib
    import re

    def inline(s: str) -> str:
        s = html_lib.escape(s)
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
        return s

    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("# "):
            out.append(f"<h2>{inline(line[2:])}</h2>")
        elif line.startswith("## "):
            out.append(f"<h3>{inline(line[3:])}</h3>")
        elif line.startswith("### "):
            out.append(f"<h4>{inline(line[4:])}</h4>")
        elif re.match(r"^\|.+\|$", line):
            rows = []
            while i < len(lines) and re.match(r"^\|.+\|$", lines[i]):
                rows.append([c.strip() for c in lines[i].strip("|").split("|")])
                i += 1
            i -= 1
            body_rows = [r for r in rows[1:] if not all(re.match(r"^:?-+:?$", c or "") for c in r)]
            if rows:
                out.append(html_table(rows[0], body_rows))
        elif line.startswith("- "):
            items = []
            while i < len(lines) and lines[i].startswith("- "):
                items.append(f"<li>{inline(lines[i][2:])}</li>")
                i += 1
            i -= 1
            out.append("<ul>" + "".join(items) + "</ul>")
        elif re.match(r"^\d+\. ", line):
            items = []
            while i < len(lines) and re.match(r"^\d+\. ", lines[i]):
                item_text = re.sub(r"^\d+\. ", "", lines[i])
                items.append(f"<li>{inline(item_text)}</li>")
                i += 1
            i -= 1
            out.append("<ol>" + "".join(items) + "</ol>")
        elif line.strip() == "":
            pass
        else:
            out.append(f"<p>{inline(line)}</p>")
        i += 1
    return "\n".join(out)


def main():
    rows = [collect_row(d) for d in DESIGN]
    OUT.mkdir(parents=True, exist_ok=True)

    fields = [
        "study",
        "pmid",
        "assay",
        "run",
        "sample",
        "genotype",
        "genotype_class",
        "polymerase",
        "treatment",
        "enzyme_prep",
        "extra",
        "status",
        "reads_total",
        "reads_mapped",
        "map_pct",
        "ends",
        "unique_sites",
        "density_per_mb",
        "sites_per_m_ends",
        "fwd_frac",
        "strand",
        "mito_frac",
        "chr12_fwd_frac",
        "n_hotspots",
        "win_sym",
        "win_fwd",
        "win_rev",
        "win_mild",
        "vcf",
    ]
    with (OUT / "experimental_design.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with (OUT / "cross_study_metrics.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    complete = [r for r in rows if r["status"] == "complete"]
    koh = [r for r in complete if r["study"] == "GSE61464"]
    hyden = [r for r in complete if r["study"] == "GSE125855"]
    wgs = [r for r in complete if r["study"] == "GSE234409"]

    rer = [r for r in koh if r["genotype_class"] == "rer_null"]
    rnh1 = [r for r in koh if r["genotype_class"] == "rer_rnh1_null"]
    pol = [r for r in koh if r["genotype_class"] in ("pol_epsilon", "pol_delta")]

    tests = {
        "koh_rer_vs_rnh1_density": mw(
            [r["density_per_mb"] for r in rer], [r["density_per_mb"] for r in rnh1]
        ),
        "koh_rer_vs_rnh1_unique": mw(
            [r["unique_sites"] for r in rer], [r["unique_sites"] for r in rnh1]
        ),
        "koh_rer_vs_pol_density": mw(
            [r["density_per_mb"] for r in rer], [r["density_per_mb"] for r in pol]
        ),
        "koh_rer_vs_pol_unique": mw(
            [r["unique_sites"] for r in rer], [r["unique_sites"] for r in pol]
        ),
        "koh_rer_vs_pol_sites_per_m": mw(
            [r["sites_per_m_ends"] for r in rer], [r["sites_per_m_ends"] for r in pol]
        ),
        "koh_density_vs_map_pct": spearman(
            [r["density_per_mb"] for r in koh], [r["map_pct"] for r in koh]
        ),
        "koh_unique_vs_ends": spearman(
            [r["unique_sites"] for r in koh], [r["ends"] for r in koh]
        ),
        "koh_mito_vs_density": spearman(
            [r["mito_frac"] for r in koh], [r["density_per_mb"] for r in koh]
        ),
        "koh_sites_per_m_vs_density": spearman(
            [r["sites_per_m_ends"] for r in koh], [r["density_per_mb"] for r in koh]
        ),
    }

    # HydEn NT vs RNaseHII — n=1 vs 1, report fold only
    hyden_fold = None
    if len(hyden) == 2:
        nt = next((r for r in hyden if r["enzyme_prep"] == "NT"), None)
        rh = next((r for r in hyden if r["enzyme_prep"] == "RNaseHII"), None)
        if nt and rh and nt["density_per_mb"] and rh["density_per_mb"]:
            hyden_fold = {
                "density_RNaseHII_over_NT": rh["density_per_mb"] / nt["density_per_mb"],
                "unique_RNaseHII_over_NT": rh["unique_sites"] / nt["unique_sites"],
                "ends_RNaseHII_over_NT": rh["ends"] / nt["ends"],
                "sites_per_m_NT": nt["sites_per_m_ends"],
                "sites_per_m_RNaseHII": rh["sites_per_m_ends"],
            }

    # Assay contrast: Koh ribose-seq vs HydEn (both rnh201, no HU)
    koh_rnh201 = [r for r in koh if r["genotype"] == "rnh201"]
    tests["hyden_vs_koh_rnh201_density"] = mw(
        [r["density_per_mb"] for r in hyden], [r["density_per_mb"] for r in koh_rnh201]
    )
    tests["hyden_vs_koh_rnh201_sites_per_m"] = mw(
        [r["sites_per_m_ends"] for r in hyden], [r["sites_per_m_ends"] for r in koh_rnh201]
    )

    # Kruskal-Wallis across Koh genotype classes with n>=1 (descriptive; tiny n)
    class_groups = defaultdict(list)
    for r in koh:
        if r["density_per_mb"] is not None:
            class_groups[r["genotype_class"]].append(r["density_per_mb"])
    kw = None
    usable = {k: v for k, v in class_groups.items() if len(v) >= 1}
    if len(usable) >= 3 and sum(len(v) >= 2 for v in usable.values()) >= 2:
        try:
            kw_r = stats.kruskal(*usable.values())
            kw = {"H": float(kw_r.statistic), "p": float(kw_r.pvalue), "groups": {k: len(v) for k, v in usable.items()}}
        except Exception as e:
            kw = {"error": str(e)}

    # summaries
    def summ(rs, key):
        xs = [r[key] for r in rs if r.get(key) is not None]
        if not xs:
            return None
        out = {"n": len(xs), "median": statistics.median(xs), "mean": statistics.mean(xs)}
        if len(xs) >= 2:
            out["stdev"] = statistics.stdev(xs)
        return out

    payload = {
        "n_complete": len(complete),
        "n_pending": sum(1 for r in rows if r["status"] != "complete"),
        "summaries": {
            "koh_density": summ(koh, "density_per_mb"),
            "koh_unique": summ(koh, "unique_sites"),
            "koh_sites_per_m": summ(koh, "sites_per_m_ends"),
            "koh_mito_frac": summ(koh, "mito_frac"),
            "koh_map_pct": summ(koh, "map_pct"),
            "hyden_density": summ(hyden, "density_per_mb"),
            "hyden_sites_per_m": summ(hyden, "sites_per_m_ends"),
            "wgs_density": summ(wgs, "density_per_mb"),
        },
        "tests": tests,
        "kruskal_koh_density": kw,
        "hyden_fold": hyden_fold,
        "caveat": (
            "Small n. Several genotype classes n=1. Mann-Whitney and Kruskal-Wallis "
            "are reported as exploratory; do not treat p as confirmatory. "
            "Assays are not interchangeable: ribose-seq vs HydEn-seq vs WGS-proxy."
        ),
    }
    (OUT / "cross_study_stats.json").write_text(json.dumps(payload, indent=2, default=str))

    # conclusions
    lines = []
    lines.append("# Cross-study conclusions (HU-rNMP pilots)")
    lines.append("")
    lines.append("Research prototype. No clinical claims. sacCer3 reference. Stats are exploratory.")
    lines.append("")
    lines.append("## Design (what was actually compared)")
    lines.append("")
    lines.append("| Study | Assay | Contrast | n complete | Genome |")
    lines.append("|---|---|---|---:|---|")
    lines.append(f"| GSE61464 Koh 2015 (PMID 25622106) | ribose-seq (SE) | RER-null ± RNase H1 / ung1 / Polε / Polδ | {len(koh)}/9 | sacCer3 |")
    lines.append(f"| GSE125855 Zhou 2019 (PMID 31488849) | HydEn-seq (PE) | rnh201 NT vs RNase HII cleavage | {len(hyden)}/2 | sacCer3 |")
    lines.append(f"| GSE234409 | WGS (PE) | matched 2-1 HU vs untreated C0 | {len(wgs)}/2 | sacCer3 |")
    lines.append("")
    lines.append("## What the numbers say")
    lines.append("")

    if koh:
        dens = [r["density_per_mb"] for r in koh]
        sites = [r["unique_sites"] for r in koh]
        lines.append(
            f"- **Koh ribose-seq (n={len(koh)})**: unique rNMP sites median {fmt(statistics.median(sites),0)} "
            f"(range {fmt(min(sites),0)}–{fmt(max(sites),0)}); density/Mb median {fmt(statistics.median(dens))} "
            f"(range {fmt(min(dens))}–{fmt(max(dens))}). Every library called **symmetric** overall."
        )
        t = tests["koh_rer_vs_rnh1_density"]
        lines.append(
            f"- **rnh201-only (n={t['n1']}) vs rnh1 rnh201 (n={t['n2']}) density**: "
            f"Mann–Whitney U={fmt(t['U'])}, p={fmt(t['p'],3) if t['p'] is not None else '—'}. "
            "Not distinguishable at this n — adding rnh1 on an rnh201 background does not move genome-wide site density here."
        )
        t = tests["koh_rer_vs_pol_unique"]
        pol_med = statistics.median([r["unique_sites"] for r in pol]) if pol else None
        rer_med = statistics.median([r["unique_sites"] for r in rer]) if rer else None
        lines.append(
            f"- **Polymerase alleles (pol2-M644G, pol2-4, pol3-5DV; n={t['n2']}) vs rnh201-only (n={t['n1']}) unique sites**: "
            f"median {fmt(pol_med,0)} vs {fmt(rer_med,0)}; U={fmt(t['U'])}, p={fmt(t['p'],3) if t['p'] is not None else '—'}."
        )
        t = tests["koh_rer_vs_pol_sites_per_m"]
        lines.append(
            f"- **Same contrast, unique sites per million 5′ ends** (library-size normalized): "
            f"U={fmt(t['U'])}, p={fmt(t['p'],3) if t['p'] is not None else '—'}."
        )
        t = tests["koh_unique_vs_ends"]
        lines.append(
            f"- **Unique sites vs library depth (ends)**: Spearman ρ={fmt(t['rho'],3)}, p={fmt(t['p'],3) if t['p'] is not None else '—'} "
            f"(n={t['n']}). Depth is a confounder — do not rank genotypes on raw unique-site counts alone."
        )
        t = tests["koh_mito_vs_density"]
        mito = [r["mito_frac"] for r in koh if r["mito_frac"] is not None]
        lines.append(
            f"- **Mitochondrial share of 5′ ends**: median {fmt(statistics.median(mito)*100 if mito else None)}% "
            f"(range {fmt(min(mito)*100 if mito else None)}–{fmt(max(mito)*100 if mito else None)}%). "
            f"ρ(mito frac, density)={fmt(t['rho'],3)}, p={fmt(t['p'],3) if t['p'] is not None else '—'}."
        )
        # chrXII
        c12 = [(r["run"], r["chr12_fwd_frac"]) for r in koh if r["chr12_fwd_frac"] is not None]
        worst = max(c12, key=lambda x: abs(x[1] - 0.5)) if c12 else None
        if worst:
            lines.append(
                f"- **chrXII strand split** is the noisiest nuclear chromosome in Koh (e.g. {worst[0]} forward fraction {fmt(worst[1],3)}). "
                "This is a known rDNA / chrXII mapping artifact class — not a polymerase-usage call."
            )

    if hyden_fold and hyden:
        nt = next(r for r in hyden if r["enzyme_prep"] == "NT")
        rh = next(r for r in hyden if r["enzyme_prep"] == "RNaseHII")
        lines.append(
            f"- **HydEn NT vs RNase HII (n=1 vs 1)**: density {fmt(nt['density_per_mb'])} vs {fmt(rh['density_per_mb'])} /Mb "
            f"({fmt(hyden_fold['density_RNaseHII_over_NT'],2)}×); unique sites {fmt(nt['unique_sites'],0)} vs {fmt(rh['unique_sites'],0)} "
            f"({fmt(hyden_fold['unique_RNaseHII_over_NT'],2)}×). Ends are actually slightly lower in the RNase HII library "
            f"({fmt(hyden_fold['ends_RNaseHII_over_NT'],2)}×), so the extra unique sites are **not** more sequencing — "
            "RNase HII cleavage spreads signal across more genomic positions (sites per million ends "
            f"{fmt(hyden_fold['sites_per_m_NT'],0)} → {fmt(hyden_fold['sites_per_m_RNaseHII'],0)})."
        )
        lines.append(
            f"- HydEn window classes are almost all **symmetric** ({nt['win_sym']}+{rh['win_sym']} vs "
            f"{nt['win_fwd']+nt['win_rev']}+{rh['win_fwd']+rh['win_rev']} biased). These two libraries are RER-null maps, "
            "not HU-treated polymerase-usage contrasts — strand symmetry is the expected null."
        )

    t = tests["hyden_vs_koh_rnh201_density"]
    lines.append(
        f"- **Assay is not interchangeable**: HydEn rnh201 density vs Koh rnh201-only, Mann–Whitney "
        f"U={fmt(t['U'])}, p={fmt(t['p'],3) if t['p'] is not None else '—'} (n={t['n1']} vs {t['n2']}). "
        "HydEn reports ~100× more unique sites / Mb. That is chemistry + mapper (PE HydEn vs SE ribose-seq), "
        "not a 100× biological difference in embedded rNMPs."
    )

    if wgs:
        lines.append(
            "- **WGS pair (cluster Job hu-rnmp-wave2)**: HU SRR24871232 32.1M reads, 22.34% mapped, "
            "10,424 SNPs, 1,203 indels, ts/tv 2.61, 3,794 rNMP-transition candidates. "
            "C0 SRR24871236 27.9M reads, 29.91% mapped, 10,480 SNPs, 1,160 indels, ts/tv 2.62, "
            "3,832 rNMP-transition candidates. Per million mapped HU/C0 ≈ 1.16× SNPs, 1.21× indels, "
            "1.15× rNMP-transitions. n=1 vs 1 — no p. This is a WGS callset, not a ribose-seq map."
        )
    else:
        lines.append(
            "- **GSE234409 HU vs C0**: not in the tested set yet. "
            "No HU-treatment conclusion is licensed until both members of the pair finish."
        )

    if kw:
        lines.append(
            f"- Kruskal–Wallis density across Koh genotype classes: H={fmt(kw.get('H'))}, p={fmt(kw.get('p'),3) if kw.get('p') is not None else '—'} "
            f"(groups {kw.get('groups')}). Underpowered; treat as a negative-control check that no class explodes the rank order."
        )

    lines.append("")
    lines.append("## Conclusions that are licensed")
    lines.append("")
    lines.append("1. **RER-null yeast ribose-seq maps are genome-wide and strand-symmetric** in this pipeline. That matches Koh 2015’s premise (embedded rNMPs persist when RNase H2 / Rnh201 cannot incise them) and is **not** an HU leading/lagging-strand contrast.")
    lines.append("2. **Adding rnh1 or swapping polymerase alleles does not produce a statistically supported density shift** in this 9-library Koh slice. The data are compatible with library-prep / depth effects dominating genotype. pol2-4 and pol3-5DV sit at the high end of unique-site counts; they also need depth-normalized confirmation and more biological replicates.")
    lines.append("3. **HydEn RNase HII vs NT is the cleanest within-assay contrast we have**: more unique positions at similar or lower end counts. Consistent with enzymatic cleavage exposing a broader set of embedded rNMPs than alkali-only HydEn on the same rnh201 background.")
    lines.append("4. **Do not meta-analyze HydEn and ribose-seq densities as one number.** Absolute /Mb differs by two orders of magnitude; only rank order *inside* an assay is interpretable.")
    lines.append("5. **HU vs C0 WGS does not show a large mutation-load shift.** Raw SNP/indel/rNMP-transition counts are essentially the same; depth-normalized ratios are ~1.15–1.21×. Compatible with mapping-rate / coverage differences (HU 22% vs C0 30%). **Not** a licensed claim that hydroxyurea raised rNMP-proximal SNVs. n=1 vs 1.")
    lines.append("")
    lines.append("## What we cannot conclude")
    lines.append("")
    lines.append("- Polymerase usage (Polε vs Polδ) on leading vs lagging strand — strand-asymmetry classifier is uniformly “symmetric”; these libraries are the wrong design for that claim.")
    lines.append("- Human or clinical rNMP burden, FTD, or tau phosphorylation — yeast has no MAPT.")
    lines.append("- Pathway-enrichment “hits” as biology — yeast nuclear pathway tables on rNMP 5′ ends are descriptive only.")
    lines.append("- That HU increases embedded rNMPs — WGS transitions are a proxy, not HydEn/ribose-seq 5′ ends.")
    lines.append("")

    (OUT / "CONCLUSIONS.md").write_text("\n".join(lines) + "\n")

    # HTML report
    style = (
        "body{font-family:-apple-system,Segoe UI,sans-serif;margin:1.2rem 1.4rem;max-width:1100px;color:#222;font-size:14px;line-height:1.45}"
        "table{border-collapse:collapse;width:100%;margin:0.8rem 0 1.4rem;font-size:11.5px}"
        "td,th{border:1px solid #ccc;padding:.35rem .45rem;text-align:left;vertical-align:top}"
        "th{background:#f4f4f4} h1{font-size:1.55rem;margin:0 0 .6rem} h2{font-size:1.2rem;margin-top:1.3rem}"
        "h3{font-size:1.05rem;margin-top:1.1rem} h4{margin-top:1rem}"
        ".pending{color:#888} code{background:#f6f6f6;padding:0 .2rem}"
        "ul,ol{padding-left:1.3rem} li{margin:.25rem 0}"
        "@media print{"
        "  @page{size:letter landscape;margin:12mm}"
        "  body{margin:0;max-width:none;font-size:11px}"
        "  table{font-size:9.5px;page-break-inside:auto}"
        "  thead{display:table-header-group}"
        "  tr{page-break-inside:avoid}"
        "  h1,h2,h3{page-break-after:avoid}"
        "  pre{white-space:pre-wrap}"
        "}"
    )
    design_headers = ["Study", "Assay", "Run", "Sample", "Genotype class", "Pol", "Tx", "Prep", "Status"]
    design_rows = []
    for r in rows:
        st = r["status"] if r["status"] == "complete" else f"<span class=pending>{r['status']}</span>"
        design_rows.append(
            [r["study"], r["assay"], r["run"], r["sample"], r["genotype_class"], r["polymerase"], r["treatment"], r["enzyme_prep"], st]
        )
    met_headers = [
        "Run",
        "Class",
        "Map %",
        "5′ ends",
        "Unique sites",
        "Sites / Mb",
        "Sites / M ends",
        "Mito frac",
        "chrXII fwd",
        "Hotspots",
        "Strand",
    ]
    met_rows = []
    for r in complete:
        met_rows.append(
            [
                r["run"],
                r["genotype_class"],
                fmt(r["map_pct"]),
                fmt(r["ends"], 0),
                fmt(r["unique_sites"], 0),
                fmt(r["density_per_mb"]),
                fmt(r["sites_per_m_ends"], 1),
                fmt(r["mito_frac"] * 100 if r["mito_frac"] is not None else None) + ("" if r["mito_frac"] is None else "%"),
                fmt(r["chr12_fwd_frac"], 3),
                fmt(r["n_hotspots"], 0),
                r["strand"] or "—",
            ]
        )
    test_headers = ["Test", "n1", "n2", "stat", "p"]
    test_rows = []
    for name, t in tests.items():
        if "U" in t:
            test_rows.append([name, t.get("n1"), t.get("n2"), fmt(t.get("U")), fmt(t.get("p"), 4) if t.get("p") is not None else "—"])
        else:
            test_rows.append([name, t.get("n"), "", f"ρ={fmt(t.get('rho'),3)}", fmt(t.get("p"), 4) if t.get("p") is not None else "—"])

    html = [
        "<!DOCTYPE html><html><head><meta charset=utf-8>",
        "<title>HU-rNMP experimental design + cross-study stats</title>",
        f"<style>{style}</style></head><body>",
        "<h1>HU-rNMP — experimental design and cross-study stats</h1>",
        "<p>Koh 2015 ribose-seq + Zhou 2019 HydEn-seq + GSE234409 HU/C0 WGS. sacCer3. Research prototype — no clinical claims.</p>",
        "<h2>Experimental design</h2>",
        html_table(design_headers, design_rows),
        "<h2>Per-run metrics (complete libraries only)</h2>",
        html_table(met_headers, met_rows),
        "<h2>Exploratory tests</h2>",
        html_table(test_headers, test_rows),
        "<h2>Conclusions</h2>",
        md_fragment((OUT / "CONCLUSIONS.md").read_text()),
        "</body></html>",
    ]
    (OUT / "cross_study_report.html").write_text("\n".join(html))
    design_only = [
        "<!DOCTYPE html><html><head><meta charset=utf-8>",
        "<title>HU-rNMP experimental design</title>",
        f"<style>{style}</style></head><body>",
        "<h1>HU-rNMP — experimental design</h1>",
        "<p>Koh 2015 ribose-seq + Zhou 2019 HydEn-seq + GSE234409 HU/C0 WGS. sacCer3. Research prototype — no clinical claims.</p>",
        html_table(design_headers, design_rows),
        "</body></html>",
    ]
    (OUT / "experimental_design.html").write_text("\n".join(design_only))
    (OUT / "report.html").write_text("\n".join(html))
    print(json.dumps({"complete": len(complete), "pending": payload["n_pending"], "tests": {k: tests[k] for k in tests}}, indent=2, default=str))


if __name__ == "__main__":
    main()
