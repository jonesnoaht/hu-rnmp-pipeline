#!/usr/bin/env python3
"""Intersect WGS T→C/A→G candidates with HydEn hotspots and transcription.

Licensed reading is a spatial co-occurrence on sacCer3. HydEn hotspots are
Zhou 2019 rnh201 maps (SRR8503261 NT, SRR8503270 RNase HII). WGS candidates
are the GSE234409 pair (SRR24871232 HU, SRR24871236 C0), n=1 vs 1. That pair
is not a licensed hydroxyurea effect. rnmp_mutations.tsv is a type/count
summary and is not read.

Coordinate contract:
  WGS pos is 1-based VCF. Hotspot bins are 0-based half-open 500 bp windows
  from (pysam_pos // 500) * 500. A VCF position overlaps a bin when
  start <= (pos - 1) < end. VCF pos 500 therefore stays in [0, 500).

Continued-stall window (all three):
  1. The same 500 bp bin is a hotspot in both HydEn libraries.
  2. At least one HU or C0 T→C/A→G candidate falls in that bin.
  3. At least one of those candidates sits in a gene body or a 500 bp
     strand-aware promoter. Intergenic-only bins are listed and not flagged.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

HYDEN_NT = "SRR8503261"
HYDEN_HII = "SRR8503270"
HU = "SRR24871232"
C0 = "SRR24871236"
WINDOW = 500
PROMOTER_UPSTREAM = 500
TRANSITIONS = {("T", "C"), ("A", "G")}
GENE_TYPES = {
    "gene",
    "ncRNA_gene",
    "rRNA_gene",
    "tRNA_gene",
    "snRNA_gene",
    "snoRNA_gene",
}

CHROM_ALIAS = {
    "i": "chrI",
    "ii": "chrII",
    "iii": "chrIII",
    "iv": "chrIV",
    "v": "chrV",
    "vi": "chrVI",
    "vii": "chrVII",
    "viii": "chrVIII",
    "ix": "chrIX",
    "x": "chrX",
    "xi": "chrXI",
    "xii": "chrXII",
    "xiii": "chrXIII",
    "xiv": "chrXIV",
    "xv": "chrXV",
    "xvi": "chrXVI",
    "m": "chrmt",
    "mt": "chrmt",
    "chrmito": "chrmt",
    "mito": "chrmt",
    "mitochondrion": "chrmt",
    "2-micron": "2micron",
    "2micron": "2micron",
}


def norm_chrom(raw: str) -> str:
    text = (raw or "").strip()
    if not text:
        return text
    if text.lower() in {"chrmt", "chrmito"}:
        return "chrmt"
    key = text.lower().replace("chromosome_", "").replace("chromosome", "").replace("chr", "")
    key = key.replace("_", "").strip()
    if key in CHROM_ALIAS:
        return CHROM_ALIAS[key]
    return text if text.startswith("chr") else f"chr{text}"


def wgs_pos_to_zero(vcf_pos: int) -> int:
    return vcf_pos - 1


def bin_start(pos0: int, window: int = WINDOW) -> int:
    if pos0 < 0:
        raise ValueError(f"0-based position must be >= 0, got {pos0}")
    return (pos0 // window) * window


def promoter_span(gene: dict, upstream: int = PROMOTER_UPSTREAM) -> tuple[int, int] | None:
    if gene.get("strand") == "+":
        return gene["start"] - upstream, gene["start"]
    if gene.get("strand") == "-":
        return gene["end"], gene["end"] + upstream
    return None


def _parse_attrs(field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in field.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def parse_sequence_regions(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("##sequence-region"):
                if line.startswith("##FASTA"):
                    break
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            lengths[norm_chrom(parts[1])] = int(parts[3])
    return lengths


def parse_genes(path: Path) -> list[dict]:
    genes = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in GENE_TYPES:
                continue
            attrs = _parse_attrs(parts[8])
            systematic = (attrs.get("gene_id") or attrs.get("ID") or "").replace("gene:", "")
            symbol = attrs.get("Name") or attrs.get("gene") or systematic
            genes.append(
                {
                    "chrom": norm_chrom(parts[0]),
                    "start": int(parts[3]) - 1,
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "systematic": systematic,
                    "symbol": symbol,
                    "biotype": attrs.get("biotype") or parts[2],
                }
            )
    return genes


def load_hotspots(path: Path) -> dict[tuple[str, int], dict]:
    found: dict[tuple[str, int], dict] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = norm_chrom(row.get("chrom") or row.get("chr") or "")
            start = int(float(row["start"]))
            end = int(float(row["end"]))
            rec = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "total": float(row.get("total") or 0),
                "fold": float(row.get("fold_over_median") or 0),
            }
            key = (chrom, start)
            prev = found.get(key)
            if prev is None or rec["total"] > prev["total"]:
                found[key] = rec
    return found


def load_wgs_transitions(path: Path, run: str) -> list[dict]:
    sites = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ref = (row.get("ref") or "").upper()
            alt = (row.get("alt") or "").upper()
            kind = row.get("type") or ""
            if kind != "rNMP_transition" or (ref, alt) not in TRANSITIONS:
                continue
            vcf_pos = int(float(row["pos"]))
            if vcf_pos < 1:
                continue
            sites.append(
                {
                    "run": run,
                    "chrom": norm_chrom(row.get("chrom") or ""),
                    "pos": vcf_pos,
                    "pos0": wgs_pos_to_zero(vcf_pos),
                    "ref": ref,
                    "alt": alt,
                }
            )
    return sites


def index_genes(genes: list[dict]) -> dict[str, list[dict]]:
    by: dict[str, list[dict]] = defaultdict(list)
    for gene in genes:
        by[gene["chrom"]].append(gene)
    for chrom in by:
        by[chrom].sort(key=lambda gene: (gene["start"], gene["end"], gene["systematic"]))
    return by


def site_transcription(pos0: int, genes: list[dict]) -> tuple[str, list[dict], list[dict]]:
    bodies = [gene for gene in genes if gene["start"] <= pos0 < gene["end"]]
    promoters = []
    for gene in genes:
        span = promoter_span(gene)
        if span and span[0] <= pos0 < span[1]:
            promoters.append(gene)
    if bodies and promoters:
        klass = "both"
    elif bodies:
        klass = "gene_body"
    elif promoters:
        klass = "promoter"
    else:
        klass = "intergenic"
    return klass, bodies, promoters


def _ids(genes: list[dict]) -> str:
    return ";".join(f"{gene['systematic']}|{gene['symbol']}" for gene in genes)


def hotspot_class(in_nt: bool, in_hii: bool) -> str:
    if in_nt and in_hii:
        return "continued"
    if in_hii:
        return "hii_only"
    if in_nt:
        return "nt_only"
    return "none"


def annotate_sites(
    sites: list[dict],
    nt: dict[tuple[str, int], dict],
    hii: dict[tuple[str, int], dict],
    genes_by_chrom: dict[str, list[dict]],
    shared_keys: set[tuple],
) -> list[dict]:
    rows = []
    for site in sites:
        start = bin_start(site["pos0"])
        end = start + WINDOW
        nt_hit = nt.get((site["chrom"], start))
        hii_hit = hii.get((site["chrom"], start))
        # A grid key can exist while the recorded end differs. Require the
        # 0-based position to fall inside the recorded half-open interval.
        in_nt = bool(nt_hit and nt_hit["start"] <= site["pos0"] < nt_hit["end"])
        in_hii = bool(hii_hit and hii_hit["start"] <= site["pos0"] < hii_hit["end"])
        klass = hotspot_class(in_nt, in_hii)
        tx, bodies, promoters = site_transcription(site["pos0"], genes_by_chrom.get(site["chrom"], []))
        key = (site["chrom"], site["pos"], site["ref"], site["alt"])
        rows.append(
            {
                **site,
                "bin_start": start,
                "bin_end": end,
                "in_nt_hotspot": int(in_nt),
                "in_hii_hotspot": int(in_hii),
                "hotspot_class": klass,
                "nt_total": nt_hit["total"] if in_nt and nt_hit else 0,
                "nt_fold": nt_hit["fold"] if in_nt and nt_hit else 0,
                "hii_total": hii_hit["total"] if in_hii and hii_hit else 0,
                "hii_fold": hii_hit["fold"] if in_hii and hii_hit else 0,
                "tx_class": tx,
                "gene_body": _ids(bodies),
                "promoter": _ids(promoters),
                "protein_coding": int(
                    any(gene["biotype"] == "protein_coding" for gene in bodies + promoters)
                ),
                "shared_site": int(key in shared_keys),
                "continued_stall": int(klass == "continued" and tx in {"gene_body", "promoter", "both"}),
            }
        )
    return rows


def end_proximal(chrom: str, start: int, end: int, chrom_lengths: dict[str, int], flank: int = 20_000) -> int:
    """Yeast chromosome-end occupancy. 20 kb is a descriptive flank, not a filter."""
    if chrom == "chrmt":
        return 1
    length = chrom_lengths.get(chrom)
    if not length:
        return 0
    return int(start < flank or end > length - flank)


def build_windows(annotated: list[dict], chrom_lengths: dict[str, int] | None = None) -> list[dict]:
    grouped: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in annotated:
        if row["hotspot_class"] == "none":
            continue
        grouped[(row["chrom"], row["bin_start"])].append(row)
    windows = []
    for (chrom, start), members in grouped.items():
        hu = [row for row in members if row["run"] == HU]
        c0 = [row for row in members if row["run"] == C0]
        shared = [row for row in members if row["shared_site"]]
        # shared is counted twice (once per arm). Unique keys:
        shared_keys = {(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in shared}
        tx_counts = defaultdict(int)
        gene_ids = []
        protein = 0
        for row in members:
            tx_counts[row["tx_class"]] += 1
            if row["gene_body"]:
                gene_ids.extend(row["gene_body"].split(";"))
            if row["promoter"]:
                gene_ids.extend(row["promoter"].split(";"))
            protein = max(protein, row["protein_coding"])
        tx_present = {row["tx_class"] for row in members}
        transcribed = tx_present & {"gene_body", "promoter", "both"}
        if transcribed == {"gene_body"}:
            tx_class = "gene_body"
        elif transcribed == {"promoter"}:
            tx_class = "promoter"
        elif len(transcribed) > 1 or "both" in transcribed:
            tx_class = "mixed" if len(transcribed) > 1 else "both"
        else:
            tx_class = "intergenic"
        klass = members[0]["hotspot_class"]
        n_hu = len(hu)
        n_c0 = len(c0)
        continued = int(klass == "continued" and bool(transcribed))
        windows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": members[0]["bin_end"],
                "hotspot_class": klass,
                "nt_total": members[0]["nt_total"],
                "nt_fold": members[0]["nt_fold"],
                "hii_total": members[0]["hii_total"],
                "hii_fold": members[0]["hii_fold"],
                "n_hu": n_hu,
                "n_c0": n_c0,
                "n_shared": len(shared_keys),
                "n_hu_only": n_hu - len({(r["pos"], r["ref"], r["alt"]) for r in hu if r["shared_site"]}),
                "n_c0_only": n_c0 - len({(r["pos"], r["ref"], r["alt"]) for r in c0 if r["shared_site"]}),
                "n_gene_body_sites": tx_counts["gene_body"] + tx_counts["both"],
                "n_promoter_sites": tx_counts["promoter"] + tx_counts["both"],
                "n_intergenic_sites": tx_counts["intergenic"],
                "tx_class": tx_class,
                "protein_coding": protein,
                "continued_stall": continued,
                "end_proximal_20kb": end_proximal(chrom, start, members[0]["bin_end"], chrom_lengths or {}),
                "genes": ";".join(sorted(set(gene for gene in gene_ids if gene))),
                "hu_over_c0": round((n_hu + 1) / (n_c0 + 1), 3),
                "license": "descriptive_cross_experiment_n1",
            }
        )
    windows.sort(key=lambda row: (-row["continued_stall"], -min(row["nt_fold"], row["hii_fold"] or 0), -row["n_hu"] - row["n_c0"], row["chrom"], row["start"]))
    return windows


def _merge_bp(spans: list[tuple[int, int]], limit: int) -> int:
    clipped = []
    for start, end in spans:
        start = max(0, start)
        end = min(limit, end)
        if end > start:
            clipped.append((start, end))
    if not clipped:
        return 0
    clipped.sort()
    cur_s, cur_e = clipped[0]
    total = 0
    for start, end in clipped[1:]:
        if start <= cur_e:
            cur_e = max(cur_e, end)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = start, end
    return total + (cur_e - cur_s)


def coverage_bp(spans: list[tuple[str, int, int]], chrom_lengths: dict[str, int]) -> int:
    by: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in spans:
        if chrom not in chrom_lengths:
            continue
        by[chrom].append((start, end))
    return sum(_merge_bp(items, chrom_lengths[chrom]) for chrom, items in by.items())


def _rate(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return numerator / denominator


def summarize(
    annotated: list[dict],
    windows: list[dict],
    genes: list[dict],
    chrom_lengths: dict[str, int],
    nt: dict[tuple[str, int], dict],
    hii: dict[tuple[str, int], dict],
) -> dict:
    nuclear = {chrom: length for chrom, length in chrom_lengths.items() if chrom != "chrmt"}
    nuclear_bp = sum(nuclear.values())

    def bin_spans(table: dict[tuple[str, int], dict]) -> list[tuple[str, int, int]]:
        return [(row["chrom"], row["start"], row["end"]) for row in table.values()]

    continued_keys = set(nt) & set(hii)
    continued_spans = [ (chrom, start, nt[(chrom, start)]["end"]) for chrom, start in continued_keys ]
    union_keys = set(nt) | set(hii)
    union_spans = []
    for key in union_keys:
        row = nt.get(key) or hii[key]
        union_spans.append((row["chrom"], row["start"], row["end"]))

    gene_spans = [(gene["chrom"], gene["start"], gene["end"]) for gene in genes]
    promoter_spans = []
    for gene in genes:
        span = promoter_span(gene)
        if span:
            promoter_spans.append((gene["chrom"], span[0], span[1]))

    def site_block(rows: list[dict]) -> dict:
        by_run = defaultdict(int)
        by_tx = defaultdict(int)
        by_class = defaultdict(int)
        for row in rows:
            by_run[row["run"]] += 1
            by_tx[row["tx_class"]] += 1
            by_class[row["hotspot_class"]] += 1
        return {
            "n": len(rows),
            "by_run": dict(by_run),
            "by_tx": dict(by_tx),
            "by_hotspot_class": dict(by_class),
            "continued_stall_sites": sum(row["continued_stall"] for row in rows),
            "shared_sites": sum(row["shared_site"] for row in rows),
        }

    hu_rows = [row for row in annotated if row["run"] == HU]
    c0_rows = [row for row in annotated if row["run"] == C0]
    continued_site_rows = [row for row in annotated if row["hotspot_class"] == "continued"]
    any_hot_rows = [row for row in annotated if row["hotspot_class"] != "none"]

    def fold(obs_n: int, obs_d: int, bg_n: int, bg_d: int) -> float | None:
        obs = _rate(obs_n, obs_d)
        bg = _rate(bg_n, bg_d)
        if obs is None or not bg:
            return None
        return round(obs / bg, 3)

    continued_bp = coverage_bp(continued_spans, nuclear)
    union_bp = coverage_bp(union_spans, nuclear)
    gene_bp = coverage_bp(gene_spans, nuclear)
    promoter_bp = coverage_bp(promoter_spans, nuclear)
    stalled = [row for row in windows if row["continued_stall"]]
    return {
        "definition": {
            "wgs": "T→C and A→G rows in rnmp_signal.tsv; not rnmp_mutations.tsv",
            "hotspots": "HydEn only (SRR8503261 NT, SRR8503270 RNase HII). Koh ribose-seq excluded.",
            "promoter_bp": PROMOTER_UPSTREAM,
            "window_bp": WINDOW,
            "coordinates": "VCF pos is 1-based; hotspot overlap uses pos-1",
            "continued_stall": "bin hot in both HydEn libraries AND a WGS transition in a gene body or 500 bp promoter",
            "hu_claim": "not licensed; n=1 vs 1, genome-wide depth-adjusted ratio previously ~1.15×",
        },
        "counts": {
            "genes": len(genes),
            "hydEn_nt_bins": len(nt),
            "hydEn_hii_bins": len(hii),
            "continued_hotspot_bins": len(continued_keys),
            "union_hotspot_bins": len(union_keys),
            "windows_with_transition": len(windows),
            "continued_stall_windows": len(stalled),
            "continued_stall_protein_coding": sum(row["protein_coding"] for row in stalled),
            "continued_stall_end_or_mt": sum(int(row.get("end_proximal_20kb") or 0) for row in stalled),
            "continued_stall_interior": sum(1 for row in stalled if not row.get("end_proximal_20kb")),
        },
        "sites": {
            "all": site_block(annotated),
            "hu": site_block(hu_rows),
            "c0": site_block(c0_rows),
            "in_any_hydEn_hotspot": site_block(any_hot_rows),
            "in_continued_hotspot": site_block(continued_site_rows),
        },
        "coverage_bp": {
            "nuclear": nuclear_bp,
            "continued_hotspot": continued_bp,
            "union_hotspot": union_bp,
            "gene_body": gene_bp,
            "promoter_500": promoter_bp,
        },
        "enrichment_vs_nuclear_bp": {
            "continued_hotspot_hu": fold(
                sum(row["hotspot_class"] == "continued" for row in hu_rows),
                len(hu_rows),
                continued_bp,
                nuclear_bp,
            ),
            "continued_hotspot_c0": fold(
                sum(row["hotspot_class"] == "continued" for row in c0_rows),
                len(c0_rows),
                continued_bp,
                nuclear_bp,
            ),
            "union_hotspot_hu": fold(
                sum(row["hotspot_class"] != "none" for row in hu_rows),
                len(hu_rows),
                union_bp,
                nuclear_bp,
            ),
            "union_hotspot_c0": fold(
                sum(row["hotspot_class"] != "none" for row in c0_rows),
                len(c0_rows),
                union_bp,
                nuclear_bp,
            ),
            "gene_body_all_wgs": fold(
                sum(row["tx_class"] in {"gene_body", "both"} for row in annotated),
                len(annotated),
                gene_bp,
                nuclear_bp,
            ),
            "promoter_all_wgs": fold(
                sum(row["tx_class"] in {"promoter", "both"} for row in annotated),
                len(annotated),
                promoter_bp,
                nuclear_bp,
            ),
        },
        "hu_vs_c0": {
            "n_hu": len(hu_rows),
            "n_c0": len(c0_rows),
            "ratio_hu_over_c0": round(len(hu_rows) / len(c0_rows), 3) if c0_rows else None,
            "note": "Counts are not depth-adjusted. Prior mapped-read ratio was ~1.15–1.21×. n=1 vs 1, no p.",
        },
    }


def write_tsv(path: Path, rows: list[dict], headers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_markdown(summary: dict, windows: list[dict]) -> str:
    counts = summary["counts"]
    sites = summary["sites"]
    enrich = summary["enrichment_vs_nuclear_bp"]
    hu = summary["hu_vs_c0"]
    stalled = [row for row in windows if row["continued_stall"]][:20]
    lines = [
        "# Stall / mutation windows vs transcription",
        "",
        "Spatial overlay on sacCer3. HydEn hotspots are Zhou 2019 *rnh201*",
        f"({HYDEN_NT} NT, {HYDEN_HII} RNase HII). WGS T→C/A→G candidates are",
        f"GSE234409 {HU} (HU) and {C0} (C0). These are different experiments.",
        "The overlay is co-occurrence, not a within-sample stall measurement.",
        "",
        "## Definition",
        "",
        "A continued-stall window is a 500 bp bin that is a hotspot in **both**",
        "HydEn libraries, contains at least one WGS T→C or A→G candidate, and",
        "places that candidate in a gene body or a 500 bp strand-aware promoter.",
        "HII-only bins are not called continued stall (NT=0 ratio class).",
        "Koh ribose-seq hotspots are not in this table. `rnmp_mutations.tsv`",
        "is a type/count summary and is not a site list.",
        "",
        "WGS positions are 1-based VCF. Overlap uses position − 1 so a call at",
        "VCF 500 stays in bin [0, 500).",
        "",
        "## Counts",
        "",
        f"- HydEn NT bins: {counts['hydEn_nt_bins']}",
        f"- HydEn HII bins: {counts['hydEn_hii_bins']}",
        f"- Continued hotspot bins (both libraries): {counts['continued_hotspot_bins']}",
        f"- Bins with a WGS transition inside a HydEn hotspot: {counts['windows_with_transition']}",
        f"- Continued-stall windows: {counts['continued_stall_windows']}",
        f"- of which touch a protein-coding feature: {counts['continued_stall_protein_coding']}",
        f"- within 20 kb of a chromosome end, or on mtDNA: {counts['continued_stall_end_or_mt']}",
        f"- interior nuclear continued-stall windows: {counts['continued_stall_interior']}",
        f"- WGS candidates: HU {hu['n_hu']} / C0 {hu['n_c0']} (raw ratio {hu['ratio_hu_over_c0']})",
        "",
        "## Transcription classes",
        "",
        "A site can be `both` when it sits in one gene's body and another gene's promoter.",
        "",
        "| Set | n | gene_body | promoter | both | intergenic |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, block in (
        ("All WGS", sites["all"]),
        ("Any HydEn hotspot", sites["in_any_hydEn_hotspot"]),
        ("Continued hotspot", sites["in_continued_hotspot"]),
    ):
        tx = block["by_tx"]
        lines.append(
            f"| {label} | {block['n']} | {tx.get('gene_body', 0)} | {tx.get('promoter', 0)} | {tx.get('both', 0)} | {tx.get('intergenic', 0)} |"
        )
    lines += [
        "",
        "## Enrichment vs nuclear length",
        "",
        "Fold is (site fraction) / (feature fraction of nuclear bp, mito excluded).",
        "Hotspot threshold is the pipeline's arbitrary ≥2× median. A fold near 1",
        "means the transitions sit in hotspots about as often as bin length predicts.",
        "",
        f"- Continued hotspot, HU: {enrich['continued_hotspot_hu']}",
        f"- Continued hotspot, C0: {enrich['continued_hotspot_c0']}",
        f"- Any HydEn hotspot, HU: {enrich['union_hotspot_hu']}",
        f"- Any HydEn hotspot, C0: {enrich['union_hotspot_c0']}",
        f"- Gene body, all WGS: {enrich['gene_body_all_wgs']}",
        f"- Promoter 500 bp, all WGS: {enrich['promoter_all_wgs']}",
        "",
        "## HU vs C0",
        "",
        hu["note"],
        "A window `hu_over_c0` above 1 is a count ratio with a pseudocount, not an effect.",
        "Shared coordinates are consistent with the same strain background, not de novo HU mutations.",
        "",
        "## Continued-stall windows (top 20 by min HydEn fold)",
        "",
        "The top of this list is chromosome-end and mtDNA occupancy. `end` means",
        "the bin is within 20 kb of a chromosome end, or on mtDNA. That flag does",
        "not remove the window from the continued-stall call.",
        "",
        "| chrom | start | end | NT fold | HII fold | HU | C0 | tx | end | genes |",
        "|---|---:|---:|---:|---:|---:|---:|---|---:|---|",
    ]
    for row in stalled:
        lines.append(
            f"| {row['chrom']} | {row['start']} | {row['end']} | {row['nt_fold']} | {row['hii_fold']} | {row['n_hu']} | {row['n_c0']} | {row['tx_class']} | {row.get('end_proximal_20kb', 0)} | {row['genes']} |"
        )
    if not stalled:
        lines.append("| — | — | — | — | — | — | — | — | — | none |")
    lines += [
        "",
        "## Not licensed",
        "",
        "- Hydroxyurea raised rNMP-proximal mutations or stalled transcription.",
        "- A gene in this table is an FTD, tau, or patient-drug locus.",
        "- HII-only bins are failed excision. Those stay in the window table as `hii_only` and are not continued-stall calls.",
        "- Rank order inside mtDNA, rDNA, or subtelomeres is copy-number / transcription occupancy, not a disease rank.",
        "",
    ]
    return "\n".join(lines)


def run(repo: Path, results: Path | None = None) -> dict:
    results = results or (repo / "results")
    gff = repo / "data/refs/sacCer3.ensembl.gff3"
    nt_path = results / "GSE125855" / HYDEN_NT / "rnmp_hotspots.tsv"
    hii_path = results / "GSE125855" / HYDEN_HII / "rnmp_hotspots.tsv"
    hu_path = results / "GSE234409" / HU / "rnmp_signal.tsv"
    c0_path = results / "GSE234409" / C0 / "rnmp_signal.tsv"
    missing = [str(path) for path in (gff, nt_path, hii_path, hu_path, c0_path) if not path.exists()]
    if missing:
        raise SystemExit("missing inputs:\n" + "\n".join(missing))

    genes = parse_genes(gff)
    chrom_lengths = parse_sequence_regions(gff)
    nt = load_hotspots(nt_path)
    hii = load_hotspots(hii_path)
    hu_sites = load_wgs_transitions(hu_path, HU)
    c0_sites = load_wgs_transitions(c0_path, C0)
    hu_keys = {(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in hu_sites}
    c0_keys = {(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in c0_sites}
    shared = hu_keys & c0_keys
    annotated = annotate_sites(hu_sites + c0_sites, nt, hii, index_genes(genes), shared)
    windows = build_windows(annotated, chrom_lengths)
    summary = summarize(annotated, windows, genes, chrom_lengths, nt, hii)
    summary["shared_unique_sites"] = len(shared)

    out = results / "_cross"
    site_headers = [
        "run", "chrom", "pos", "pos0", "ref", "alt", "bin_start", "bin_end",
        "in_nt_hotspot", "in_hii_hotspot", "hotspot_class",
        "nt_total", "nt_fold", "hii_total", "hii_fold",
        "tx_class", "gene_body", "promoter", "protein_coding",
        "shared_site", "continued_stall",
    ]
    window_headers = [
        "chrom", "start", "end", "hotspot_class", "nt_total", "nt_fold",
        "hii_total", "hii_fold", "n_hu", "n_c0", "n_shared", "n_hu_only",
        "n_c0_only", "n_gene_body_sites", "n_promoter_sites", "n_intergenic_sites",
        "tx_class", "protein_coding", "continued_stall", "end_proximal_20kb",
        "genes", "hu_over_c0", "license",
    ]
    write_tsv(out / "stall_mutation_sites.tsv", annotated, site_headers)
    write_tsv(out / "stall_mutation_windows.tsv", windows, window_headers)
    stalled = [row for row in windows if row["continued_stall"]]
    write_tsv(out / "continued_stall_windows.tsv", stalled, window_headers)
    (out / "stall_mutation_windows.json").write_text(json.dumps(summary, indent=2) + "\n")
    text = render_markdown(summary, windows)
    (out / "STALL_MUTATION_WINDOWS.md").write_text(text)
    docs = repo / "docs/stall-mutation-windows.md"
    docs.parent.mkdir(parents=True, exist_ok=True)
    docs.write_text(text)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("/home/noahtjones/hu-rnmp-pipeline"))
    parser.add_argument("--results", type=Path, default=None)
    args = parser.parse_args()
    summary = run(args.repo, args.results)
    print(json.dumps(summary["counts"]))
    print(json.dumps(summary["enrichment_vs_nuclear_bp"]))


if __name__ == "__main__":
    main()
