#!/usr/bin/env python3
"""Rank sacCer3 genes for rNMP incorporation and failed excision.

Plan task t-jj-gene-excision: intersect the SGD GFF with HydEn hotspots and
signal, then rank genes. The licensed failed-excision contrast is HydEn
RNase HII / NT on the same rnh201 background (SRR8503270 / SRR8503261).
That ratio is HII end-counts divided by NT end-counts. It is not a
pseudocount, and an NT=0 gene is not a ratio.

Koh ribose-seq (rnh201 vs rnh1 rnh201) is not read. That contrast is
underpowered and is a different assay.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import urllib.request
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

HYDEN_NT = "SRR8503261"
HYDEN_HII = "SRR8503270"
GFF_URLS = (
    "https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz",
    "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz",
)
GENE_TYPES = {
    "gene",
    "pseudogene",
    "ncRNA_gene",
    "tRNA_gene",
    "snoRNA_gene",
    "snRNA_gene",
    "rRNA_gene",
    "telomerase_RNA_gene",
    "transposable_element_gene",
}
EMPTY_INDEX: tuple[list[dict], list[int]] = ([], [])
# Both libraries must clear these floors before a ratio is a call.
# Below the floor the gene is sparse, not enriched.
MIN_ENDS = 100
MIN_SITES = 10
ENRICHED = 2.0
DEPLETED = 0.5
END_PROXIMAL_BP = 20_000

CHROM_ALIAS = {
    "chrmt": "chrM",
    "chrmito": "chrM",
    "mito": "chrM",
    "mitochondrion": "chrM",
    "mt": "chrM",
}


def norm_chrom(name: str) -> str:
    """Map SGD chrmt onto the pipeline's sacCer3 chrM. Other names stay."""
    key = (name or "").strip()
    return CHROM_ALIAS.get(key.lower(), key)


def parse_attrs(field: str) -> dict[str, str]:
    out = {}
    for kv in field.split(";"):
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        out[k] = unquote(v)
    return out


def parse_gff_genes(path: Path) -> tuple[list[dict], dict[str, int]]:
    """Return gene features and chromosome lengths on pipeline chrom names.

    SGD gene spans (R64-4-1 weekly GFF) were extended in February 2024 to the
    longest experimental transcript. Coordinates stay 0-based half-open.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    genes = []
    lengths = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom = norm_chrom(parts[0])
            kind = parts[2]
            if kind == "chromosome":
                lengths[chrom] = int(parts[4])
                continue
            if kind not in GENE_TYPES:
                continue
            attrs = parse_attrs(parts[8])
            systematic = attrs.get("ID") or attrs.get("Name") or ""
            symbol = attrs.get("gene") or attrs.get("Name") or systematic
            note = attrs.get("Note") or ""
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            if end0 <= start0 or not systematic:
                continue
            genes.append({
                "chrom": chrom,
                "start": start0,
                "end": end0,
                "strand": parts[6],
                "systematic": systematic,
                "symbol": symbol,
                "biotype": kind,
                "dubious": note.lower().startswith("dubious"),
                "te": kind == "transposable_element_gene",
                "rrna": kind == "rRNA_gene" or symbol.startswith("RDN"),
                "mitochondrial": chrom == "chrM",
            })
    return genes, lengths


def end_proximal(gene: dict, lengths: dict[str, int]) -> bool:
    length = lengths.get(gene["chrom"])
    if not length:
        return False
    if gene["chrom"] == "chrM":
        return False
    return gene["start"] < END_PROXIMAL_BP or gene["end"] > length - END_PROXIMAL_BP


def _by_chrom(genes: list[dict]) -> dict[str, tuple[list[dict], list[int]]]:
    by = defaultdict(list)
    for gene in genes:
        by[gene["chrom"]].append(gene)
    indexed = {}
    for chrom, chrom_genes in by.items():
        chrom_genes.sort(key=lambda g: (g["start"], g["end"], g["systematic"]))
        indexed[chrom] = (chrom_genes, [g["start"] for g in chrom_genes])
    return indexed


def genes_overlapping(indexed: tuple[list[dict], list[int]], start: int, end: int) -> list[dict]:
    """Genes whose [start, end) overlaps [start, end). Index is start-sorted."""
    sorted_genes, starts = indexed
    if not sorted_genes:
        return []
    # Every overlap starts at or before `end`. Walk back only while a long
    # gene could still cover `start` (yeast gene spans are far under 100 kb).
    right = bisect.bisect_left(starts, end)
    hits = []
    for i in range(right - 1, -1, -1):
        gene = sorted_genes[i]
        if gene["end"] <= start and gene["start"] < start - 100_000:
            break
        if gene["end"] > start and gene["start"] < end:
            hits.append(gene)
    hits.reverse()
    return hits


def hotspot_shares(window_start: int, window_end: int, total: float, genes: list[dict]) -> list[tuple[dict, float]]:
    """Assign window mass by overlap length / window length. No equal split."""
    width = window_end - window_start
    if width <= 0 or total == 0:
        return []
    shares = []
    for gene in genes:
        overlap = min(window_end, gene["end"]) - max(window_start, gene["start"])
        if overlap > 0:
            shares.append((gene, total * overlap / width))
    return shares


def classify_excision(nt_ends: float, hii_ends: float, nt_sites: int, hii_sites: int) -> tuple[str, float | None]:
    """Return (class, HII/NT unique-position ratio or None).

    Remaining rNMPs are positions, not read depth. The ratio is HII sites / NT
    sites. End counts are a depth check and do not set the class. NT sites = 0
    is HII_only, not a ratio. Pseudocounts are not applied.
    """
    if nt_sites <= 0 and hii_sites <= 0 and nt_ends <= 0 and hii_ends <= 0:
        return "absent", None
    if nt_sites <= 0 or nt_ends <= 0:
        return "HII_only", None
    if hii_sites <= 0 or hii_ends <= 0:
        return "NT_only", None
    ratio = hii_sites / nt_sites
    supported = (
        nt_ends >= MIN_ENDS
        and hii_ends >= MIN_ENDS
        and nt_sites >= MIN_SITES
        and hii_sites >= MIN_SITES
    )
    if not supported:
        return "sparse", ratio
    if ratio >= ENRICHED:
        return "HII_enriched", ratio
    if ratio <= DEPLETED:
        return "NT_enriched", ratio
    return "shared", ratio


def load_hotspots(path: Path) -> list[dict]:
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({
                "chrom": norm_chrom(row["chrom"]),
                "start": int(float(row["start"])),
                "end": int(float(row["end"])),
                "total": float(row["total"]),
            })
    return rows


def accumulate(genes: list[dict], hotspots_by_run: dict[str, list[dict]], signal_paths: dict[str, Path]) -> dict:
    indexed = _by_chrom(genes)
    per = {
        g["systematic"]: {
            "nt_ends": 0.0,
            "hii_ends": 0.0,
            "nt_sites": 0,
            "hii_sites": 0,
            "nt_hotspot": 0.0,
            "hii_hotspot": 0.0,
        }
        for g in genes
    }
    totals = {
        "nt_ends": 0.0,
        "hii_ends": 0.0,
        "nt_ends_intergenic": 0.0,
        "hii_ends_intergenic": 0.0,
        "nt_hotspot": 0.0,
        "hii_hotspot": 0.0,
        "nt_hotspot_intergenic": 0.0,
        "hii_hotspot_intergenic": 0.0,
    }
    run_field = {HYDEN_NT: "nt", HYDEN_HII: "hii"}
    for run, windows in hotspots_by_run.items():
        prefix = run_field[run]
        for window in windows:
            totals[f"{prefix}_hotspot"] += window["total"]
            hits = genes_overlapping(indexed.get(window["chrom"], EMPTY_INDEX), window["start"], window["end"])
            shares = hotspot_shares(window["start"], window["end"], window["total"], hits)
            assigned = 0.0
            for gene, share in shares:
                per[gene["systematic"]][f"{prefix}_hotspot"] += share
                assigned += share
            totals[f"{prefix}_hotspot_intergenic"] += window["total"] - assigned

    for run, path in signal_paths.items():
        prefix = run_field[run]
        seen = {g["systematic"]: set() for g in genes}
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                chrom = norm_chrom(row["chrom"])
                pos = int(float(row["pos"]))
                count = float(row["count"])
                totals[f"{prefix}_ends"] += count
                hits = genes_overlapping(indexed.get(chrom, EMPTY_INDEX), pos, pos + 1)
                if not hits:
                    totals[f"{prefix}_ends_intergenic"] += count
                    continue
                for gene in hits:
                    key = gene["systematic"]
                    per[key][f"{prefix}_ends"] += count
                    if pos not in seen[key]:
                        seen[key].add(pos)
                        per[key][f"{prefix}_sites"] += 1
    return {"per": per, "totals": totals}


def gene_rows(genes: list[dict], lengths: dict[str, int], tallies: dict) -> list[dict]:
    rows = []
    for gene in genes:
        stats = tallies[gene["systematic"]]
        length = gene["end"] - gene["start"]
        kb = length / 1000.0
        klass, ratio = classify_excision(
            stats["nt_ends"], stats["hii_ends"], stats["nt_sites"], stats["hii_sites"]
        )
        hot_ratio = None
        if stats["nt_hotspot"] > 0 and stats["hii_hotspot"] > 0:
            hot_ratio = stats["hii_hotspot"] / stats["nt_hotspot"]
        end_ratio = None
        if stats["nt_ends"] > 0 and stats["hii_ends"] > 0:
            end_ratio = stats["hii_ends"] / stats["nt_ends"]
        rows.append({
            "systematic": gene["systematic"],
            "symbol": gene["symbol"],
            "biotype": gene["biotype"],
            "chrom": gene["chrom"],
            "start0": gene["start"],
            "end0": gene["end"],
            "length_bp": length,
            "dubious": int(gene["dubious"]),
            "te": int(gene["te"]),
            "mitochondrial": int(gene["mitochondrial"]),
            "rrna": int(gene["rrna"]),
            "end_proximal": int(end_proximal(gene, lengths)),
            "nt_ends": round(stats["nt_ends"], 2),
            "hii_ends": round(stats["hii_ends"], 2),
            "nt_sites": stats["nt_sites"],
            "hii_sites": stats["hii_sites"],
            "nt_ends_per_kb": round(stats["nt_ends"] / kb, 4) if kb else 0,
            "hii_ends_per_kb": round(stats["hii_ends"] / kb, 4) if kb else 0,
            "nt_sites_per_kb": round(stats["nt_sites"] / kb, 4) if kb else 0,
            "hii_sites_per_kb": round(stats["hii_sites"] / kb, 4) if kb else 0,
            "hii_sites_over_nt": "" if ratio is None else round(ratio, 4),
            "hii_ends_over_nt": "" if end_ratio is None else round(end_ratio, 4),
            "nt_hotspot": round(stats["nt_hotspot"], 2),
            "hii_hotspot": round(stats["hii_hotspot"], 2),
            "hotspot_hii_over_nt": "" if hot_ratio is None else round(hot_ratio, 4),
            "excision_class": klass,
        })
    return rows


def nuclear_interior(row: dict) -> bool:
    return (
        row["biotype"] == "gene"
        and not row["dubious"]
        and not row["te"]
        and not row["mitochondrial"]
        and not row["rrna"]
        and not row["end_proximal"]
    )


def write_tsv(path: Path, rows: list[dict], headers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _fmt(n: float) -> str:
    if isinstance(n, float) and not n.is_integer():
        return f"{n:.4g}"
    return f"{n:,.0f}" if isinstance(n, (int, float)) else str(n)


def render_markdown(rows: list[dict], totals: dict, gff_path: Path, n_genes: int) -> str:
    by_class = defaultdict(int)
    for row in rows:
        by_class[row["excision_class"]] += 1
    both = [r for r in rows if r["hii_sites_over_nt"] != ""]
    supported = [r for r in rows if r["excision_class"] in {"HII_enriched", "NT_enriched", "shared"}]
    nuclear_supported = [r for r in supported if nuclear_interior(r)]
    enriched = [r for r in nuclear_supported if r["excision_class"] == "HII_enriched"]
    enriched.sort(key=lambda r: (-float(r["hii_sites_over_nt"]), -r["hii_sites"], r["systematic"]))
    depleted = [r for r in nuclear_supported if r["excision_class"] == "NT_enriched"]
    depleted.sort(key=lambda r: (float(r["hii_sites_over_nt"]), r["systematic"]))

    by_sites = [r for r in rows if nuclear_interior(r)]
    by_sites.sort(key=lambda r: (-r["hii_sites_per_kb"], -r["hii_sites"], r["systematic"]))
    by_ends = list(by_sites)
    by_ends.sort(key=lambda r: (-r["hii_ends_per_kb"], -r["hii_ends"], r["systematic"]))

    def share(predicate) -> tuple[int, float]:
        chosen = [r for r in rows if predicate(r)]
        return len(chosen), sum(r["hii_ends"] for r in chosen)

    compartments = [
        ("Mitochondrial", lambda r: r["mitochondrial"]),
        ("rRNA / RDN", lambda r: r["rrna"] and not r["mitochondrial"]),
        ("Transposable element", lambda r: r["te"]),
        ("Dubious ORF", lambda r: r["dubious"] and not r["mitochondrial"] and not r["te"]),
        ("End-proximal nuclear ORF", lambda r: r["end_proximal"] and r["biotype"] == "gene" and not r["dubious"]),
        ("Nuclear interior ORF", nuclear_interior),
    ]
    hii_total = totals["hii_ends"] or 1.0
    end_ratio = totals["hii_ends"] / totals["nt_ends"] if totals["nt_ends"] else float("nan")
    lines = [
        "# Preferential rNMP incorporation and failed excision by gene",
        "",
        "SGD GFF gene spans intersected with Zhou 2019 HydEn-seq on sacCer3.",
        f"Libraries: {HYDEN_NT} (rnh201 NT, no RNase HII pre-cleavage) and",
        f"{HYDEN_HII} (rnh201 RNase HII). Same rnh201 background. Koh ribose-seq",
        "is not in this table: rnh201 vs rnh1+rnh201 is underpowered and is a",
        "different assay.",
        "",
        "## What the ratio is",
        "",
        "RNase HII / NT is the within-assay contrast for rNMPs that remained",
        "embedded in rnh201 cells. The HII library records ends after in vitro",
        "RNase HII cleavage. The NT library records alkali / hydrolytic ends",
        "without that cleavage. These are separate libraries, not a before/after",
        "subtraction of one DNA sample. Pipeline coordinates already apply the",
        "1-nt Alk vs RHII offset (Ribose-Map coordinate.sh); this intersection",
        "does not shift them again.",
        "",
        "The ranked ratio is unique HII positions / unique NT positions on the",
        "SGD gene span. Read-end counts are a depth column. Genome-wide, HII",
        f"has fewer ends (ratio {end_ratio:.4f}) and more distinct positions,",
        "so an end-count ratio mostly tracks library depth. A call requires at",
        f"least {MIN_ENDS} ends and {MIN_SITES} distinct positions on each side,",
        f"and HII/NT ≥ {ENRICHED} on positions. NT sites = 0 is class `HII_only`",
        "and has no ratio. The pseudocount (HII+1)/(NT+1) is not used.",
        "",
        f"Gene models: `{gff_path.name}` ({n_genes} gene-like features). SGD",
        "assembly header R64-4-1; chromosome lengths match sacCer3 / R64-1-1,",
        "including chrM = chrmt = 85,779 bp. Gene spans are the February 2024",
        "transcript-extended parent features, not CDS-only intervals. SGD",
        "`chrmt` is matched to pipeline `chrM`.",
        "",
        "## Library totals",
        "",
        f"- NT ends: {_fmt(totals['nt_ends'])}; of which no SGD gene hit: {_fmt(totals['nt_ends_intergenic'])}",
        f"- HII ends: {_fmt(totals['hii_ends'])}; of which no SGD gene hit: {_fmt(totals['hii_ends_intergenic'])}",
        f"- Genome-wide HII/NT end ratio: {end_ratio:.4f}",
        f"- NT hotspot-window mass: {_fmt(totals['nt_hotspot'])}",
        f"- HII hotspot-window mass: {_fmt(totals['hii_hotspot'])}",
        "",
        "Gene end sums double-count a position that falls in two gene spans.",
        "Intergenic is the complement: ends whose position hits no gene.",
        "",
        "## Where HII ends sit",
        "",
        "| Compartment | Genes | HII ends | Share of library |",
        "|---|---:|---:|---:|",
    ]
    for label, pred in compartments:
        n, ends = share(pred)
        lines.append(f"| {label} | {n} | {_fmt(ends)} | {100 * ends / hii_total:.1f}% |")
    lines += [
        "",
        "Mitochondrial genes, rDNA, and subtelomeres can dominate raw end totals.",
        "Incorporation rank below is nuclear interior ORFs (biotype `gene`,",
        "not dubious, not within 20 kb of a chromosome end).",
        "",
        "## Preferential incorporation (nuclear interior ORFs)",
        "",
        "Primary rank: distinct HydEn HII positions per kb. That is where embedded",
        "rNMPs are spread, not where one site was sequenced deeply. A second table",
        "ranks HII ends per kb so pileups stay visible. Neither table is a human",
        "disease rank.",
        "",
        "### Distinct positions per kb (top 25)",
        "",
        "| Symbol | Systematic | Length | HII sites/kb | HII sites | NT sites | Site HII/NT | End HII/NT |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in by_sites[:25]:
        lines.append(
            f"| {row['symbol']} | {row['systematic']} | {row['length_bp']} | {row['hii_sites_per_kb']} | "
            f"{row['hii_sites']} | {row['nt_sites']} | {row['hii_sites_over_nt'] or '—'} | "
            f"{row['hii_ends_over_nt'] or '—'} |"
        )
    lines += [
        "",
        "### End depth per kb (top 15)",
        "",
        "High ends/kb with a site ratio near 1 is a pileup, not more incorporation sites.",
        "",
        "| Symbol | Systematic | HII ends/kb | HII ends | End HII/NT | Site HII/NT |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in by_ends[:15]:
        lines.append(
            f"| {row['symbol']} | {row['systematic']} | {row['hii_ends_per_kb']} | "
            f"{_fmt(row['hii_ends'])} | {row['hii_ends_over_nt'] or '—'} | {row['hii_sites_over_nt'] or '—'} |"
        )
    lines += [
        "",
        "## Failed excision classes",
        "",
        "Class uses the unique-position ratio. End HII/NT is stored and is not the call.",
        "",
        "| Class | Genes | Meaning |",
        "|---|---:|---|",
        f"| HII_enriched | {by_class['HII_enriched']} | both sides clear the floor and site HII/NT ≥ {ENRICHED} |",
        f"| shared | {by_class['shared']} | both sides clear the floor and {DEPLETED} < site ratio < {ENRICHED} |",
        f"| NT_enriched | {by_class['NT_enriched']} | both sides clear the floor and site HII/NT ≤ {DEPLETED} |",
        f"| sparse | {by_class['sparse']} | both sides nonzero, floor not met; ratio stored, not called |",
        f"| HII_only | {by_class['HII_only']} | NT sites or NT ends = 0; not a ratio |",
        f"| NT_only | {by_class['NT_only']} | HII sites or HII ends = 0; not a ratio |",
        f"| absent | {by_class['absent']} | no HydEn ends on the gene span |",
        "",
        f"Both sides nonzero: {len(both)} genes. Supported calls: {len(supported)}.",
        f"Nuclear interior supported calls: {len(nuclear_supported)}, of which "
        f"{len(enriched)} are HII_enriched.",
        "",
        "### Nuclear HII-enriched (top 25 by site ratio)",
        "",
        "More distinct positions after RNase HII cleavage than in the alkali-only",
        "library, on a nuclear interior ORF. In rnh201 cells that is the",
        "within-assay reading of embedded rNMPs the alkali library under-called.",
        "It is not an excision-rate measurement.",
        "",
        "| Symbol | Systematic | Length | Site HII/NT | HII sites | NT sites | End HII/NT |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in enriched[:25]:
        lines.append(
            f"| {row['symbol']} | {row['systematic']} | {row['length_bp']} | {row['hii_sites_over_nt']} | "
            f"{row['hii_sites']} | {row['nt_sites']} | {row['hii_ends_over_nt'] or '—'} |"
        )
    if not enriched:
        lines.append("| — | — | — | — | — | — | — |")
    lines += [
        "",
        "### Nuclear NT-enriched on positions",
        "",
        "Site HII/NT ≤ 0.5 would mean the HII library lost positions relative to",
        "alkali. That class is empty or short here because the HII library",
        "expands unique positions genome-wide. End-count NT bias (NUP49, LYS1,",
        "MEP2) is a depth pileup: those genes keep a site ratio near 1. Do not",
        "read either pattern as successful excision.",
        "",
    ]
    if depleted:
        lines += [
            "| Symbol | Systematic | Site HII/NT | HII sites | NT sites | End HII/NT |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for row in depleted[:15]:
            lines.append(
                f"| {row['symbol']} | {row['systematic']} | {row['hii_sites_over_nt']} | "
                f"{row['hii_sites']} | {row['nt_sites']} | {row['hii_ends_over_nt'] or '—'} |"
            )
    else:
        lines.append("No nuclear interior ORF clears the floor with site HII/NT ≤ 0.5.")
    lines += [
        "",
        "## What this does not say",
        "",
        "- No FTD, MAPT, or tau conclusion. Yeast has no MAPT.",
        "- No hydroxyurea effect. HU vs C0 is a different experiment and is not ranked here.",
        "- No Koh genotype contrast. Those libraries are not inputs.",
        "- HII_only genes, and any NT=0 pseudocount ratio, are not failed-excision calls.",
        "- An end-count ratio above 2 on mtDNA intron ORFs is depth, not the position rank.",
        "",
        "Tables: `docs/gene-excision/gene_rnmp_sgd.tsv`,",
        "`docs/gene-excision/incorporation_rank.tsv`,",
        "`docs/gene-excision/failed_excision_supported.tsv`.",
        "",
    ]
    return "\n".join(lines)


def ensure_gff(path: Path) -> Path:
    if path.exists() and path.stat().st_size > 100_000:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for url in GFF_URLS:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "hu-rnmp-pipeline/gene-excision"})
            with urllib.request.urlopen(req, timeout=180) as response, path.open("wb") as out:
                out.write(response.read())
            if path.stat().st_size > 100_000:
                return path
        except Exception as exc:  # noqa: BLE001 — try the archive mirror next
            last_error = exc
    raise RuntimeError(f"SGD GFF download failed: {last_error}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--results", type=Path, default=None)
    parser.add_argument("--gff", type=Path, default=None)
    args = parser.parse_args()
    repo = args.repo
    results = args.results or (repo / "results")
    gff = ensure_gff(args.gff or (repo / "data/refs/sacCer3.sgd.gff.gz"))
    genes, lengths = parse_gff_genes(gff)
    print(f"genes={len(genes)} chroms={len(lengths)} gff={gff}")
    hotspots = {
        HYDEN_NT: load_hotspots(results / "GSE125855" / HYDEN_NT / "rnmp_hotspots.tsv"),
        HYDEN_HII: load_hotspots(results / "GSE125855" / HYDEN_HII / "rnmp_hotspots.tsv"),
    }
    signals = {
        HYDEN_NT: results / "GSE125855" / HYDEN_NT / "rnmp_signal.tsv",
        HYDEN_HII: results / "GSE125855" / HYDEN_HII / "rnmp_signal.tsv",
    }
    tallies = accumulate(genes, hotspots, signals)
    rows = gene_rows(genes, lengths, tallies["per"])
    headers = list(rows[0].keys()) if rows else []
    out = repo / "docs" / "gene-excision"
    write_tsv(out / "gene_rnmp_sgd.tsv", rows, headers)
    interior = [r for r in rows if nuclear_interior(r)]
    interior.sort(key=lambda r: (-r["hii_sites_per_kb"], -r["hii_sites"], r["systematic"]))
    for i, row in enumerate(interior, start=1):
        row["incorporation_rank"] = i
    write_tsv(
        out / "incorporation_rank.tsv",
        interior,
        ["incorporation_rank", "symbol", "systematic", "chrom", "length_bp",
         "hii_sites_per_kb", "hii_sites", "nt_sites", "hii_sites_over_nt",
         "hii_ends_per_kb", "hii_ends_over_nt", "excision_class"],
    )
    supported = [r for r in rows if r["excision_class"] in {"HII_enriched", "NT_enriched", "shared"}]
    supported.sort(key=lambda r: (-float(r["hii_sites_over_nt"]), -r["hii_sites"], r["systematic"]))
    write_tsv(
        out / "failed_excision_supported.tsv",
        supported,
        ["symbol", "systematic", "biotype", "chrom", "length_bp", "hii_sites_over_nt",
         "hii_sites", "nt_sites", "hii_ends_over_nt", "hii_ends", "nt_ends",
         "excision_class", "dubious", "mitochondrial", "end_proximal"],
    )
    text = render_markdown(rows, tallies["totals"], gff, len(genes))
    (out / "GENE_EXCISION.md").write_text(text)
    (repo / "docs" / "gene-excision.md").write_text(text)
    print(
        f"wrote {out} supported={len(supported)} "
        f"HII_enriched={sum(r['excision_class']=='HII_enriched' for r in rows)} "
        f"HII_only={sum(r['excision_class']=='HII_only' for r in rows)}"
    )


if __name__ == "__main__":
    main()
