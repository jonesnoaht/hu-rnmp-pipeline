#!/usr/bin/env python3
"""Summarize rNMP calls.

Window classes are descriptive count ratios. v1 does not assign replicative
polymerases. The WGS headline metric is the 2–5 bp short-tandem-repeat
deletion count. Raw T>C / A>G counts stay under ``qc`` and are not a signal.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def _table(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line in handle:
            if not line.strip():
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
        return header, rows


POLYMERASE_NOTE = (
    "v1 does not assign replicative polymerases. Window classes are descriptive "
    "ratios of rNMP counts on the two strands. Replication origins are not modeled."
)
V2_POLYMERASE_NOTE = (
    "Future work / v2: polymerase assignment only for polymerase-variant strains, "
    "relative to replication origins and OK-seq (Daigaku et al. 2015, PMID 25664722)."
)
REPEAT_DELETION_NOTE = (
    "Headline metric: deletions of 2–5 bp that sit in a short tandem repeat "
    "(unit length 2–5 bp, at least two copies; homopolymer units excluded). "
    "This is the Top1-dependent mutation signature of rNMPs left in yeast DNA "
    "(Nick McElhinny et al. 2010, PMID 20729855; Kim et al. 2011, PMID 21700875)."
)


def _qc_counts(raw: dict | None) -> dict:
    raw = raw or {}
    tc = int(raw.get("T>C", 0) or 0)
    ag = int(raw.get("A>G", 0) or 0)
    return {"T>C": tc, "A>G": ag, "total": tc + ag}


def genome_bp_from_fai(path: str) -> int:
    total = 0
    with open(path) as handle:
        for line in handle:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1].isdigit():
                total += int(parts[1])
    return total


def scan_direct(
    signal_tsv,
    asymmetry_tsv,
    output_tsv,
    signature_json,
    assay,
    technique,
    genome_bp,
    stats_json,
    accession=None,
    genotype=None,
):
    _header, signal_rows = _table(signal_tsv)
    strand_counts = {"+": 0, "-": 0}
    base_counts: Counter[str] = Counter()
    chrom_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"+": 0, "-": 0})
    for row in signal_rows:
        try:
            count = int(row["count"])
        except (KeyError, ValueError):
            continue
        strand = row.get("strand", "")
        if strand in strand_counts:
            strand_counts[strand] += count
            chrom_counts[row.get("chrom", "")][strand] += count
        base = row.get("base", "NA")
        base_counts[base] += count

    _header, asym_rows = _table(asymmetry_tsv)
    class_counts: Counter[str] = Counter()
    biases = []
    for row in asym_rows:
        if "asymmetry" not in row or "strand_bias" not in row:
            continue
        class_counts[row["asymmetry"]] += 1
        biases.append(float(row["strand_bias"]))

    total_ends = strand_counts["+"] + strand_counts["-"]
    fwd_frac = strand_counts["+"] / max(total_ends, 1)
    avg_bias = sum(biases) / len(biases) if biases else None
    calls_per_mb = round(total_ends / (genome_bp / 1e6), 4) if genome_bp else None

    stats = {}
    if stats_json and Path(stats_json).exists():
        stats = json.loads(Path(stats_json).read_text())

    with open(output_tsv, "w") as handle:
        handle.write("metric\tvalue\n")
        lines = [
            f"total_5prime_ends\t{total_ends}\n",
            f"forward_ends\t{strand_counts['+']}\n",
            f"reverse_ends\t{strand_counts['-']}\n",
            f"forward_fraction\t{fwd_frac:.4f}\n",
            f"unique_rnmp_positions\t{len(signal_rows)}\n",
            f"assay\t{assay}\n",
            f"technique\t{technique or ''}\n",
        ]
        if avg_bias is not None:
            lines.append(f"avg_strand_bias\t{avg_bias:.4f}\n")
        lines.extend(f"window_{name}\t{count}\n" for name, count in sorted(class_counts.items()))
        lines.extend(f"base_{base}\t{count}\n" for base, count in sorted(base_counts.items()))
        handle.writelines(lines)

    per_chrom = {
        chrom: {"forward": strands["+"], "reverse": strands["-"]}
        for chrom, strands in sorted(chrom_counts.items())
    }
    signature = {
        "assay_type": assay,
        "technique": technique,
        "total_rnmp_5prime_ends": total_ends,
        "forward_ends": strand_counts["+"],
        "reverse_ends": strand_counts["-"],
        "forward_fraction": round(fwd_frac, 4),
        "unique_rnmp_positions": len(signal_rows),
        "base_counts": dict(sorted(base_counts.items())),
        "avg_strand_bias": None if avg_bias is None else round(avg_bias, 4),
        "window_classes": dict(class_counts),
        "per_chromosome": per_chrom,
        "rnmp_calls_per_mb": calls_per_mb,
        "rnmp_calls_per_mb_note": (
            "Raw called rNMP counts divided by reference length. Not normalized "
            "to library size or restriction-site ends. Null when no .fai length was given. "
            "The previous hard-coded unique_positions/12.1 value was removed."
        ),
        "polymerase_assignment": "experimental_not_assigned",
        "polymerase_note": POLYMERASE_NOTE,
        "v2_note": V2_POLYMERASE_NOTE,
        "accession": accession or "",
        "genotype": genotype or "",
        "call_stats": stats,
    }
    Path(signature_json).write_text(json.dumps(signature, indent=2) + "\n")
    print(
        f"{technique or assay}: {total_ends} rNMP counts, "
        f"{len(signal_rows)} positions, bases {dict(base_counts)} → {output_tsv}",
        file=sys.stderr,
    )
    return signature


def scan_wgs(signal_tsv, strand_tsv, output_tsv, signature_json, stats_json=None, accession=None, genotype=None):
    _header, signal_rows = _table(signal_tsv)
    deletions = []
    transition_rows: Counter[str] = Counter()
    for row in signal_rows:
        kind = row.get("type", "")
        if kind == "repeat_deletion":
            try:
                length = int(row.get("del_len") or 0)
                unit = int(row.get("repeat_unit") or 0)
            except ValueError:
                continue
            deletions.append({"del_len": length, "repeat_unit": unit})
        elif kind == "rNMP_transition":
            key = f"{row.get('ref', '')}>{row.get('alt', '')}"
            if key in {"T>C", "A>G"}:
                transition_rows[key] += 1

    deletion_sizes = Counter(row["del_len"] for row in deletions)
    unit_sizes = Counter(row["repeat_unit"] for row in deletions)

    qc_raw = None
    if stats_json and Path(stats_json).exists():
        stats = json.loads(Path(stats_json).read_text())
        qc_raw = (stats.get("qc") or {}).get("raw_tc_ag_counts")
    if qc_raw is None:
        qc_raw = {"T>C": transition_rows["T>C"], "A>G": transition_rows["A>G"]}
    qc = {
        "raw_tc_ag_counts": _qc_counts(qc_raw),
        "note": (
            "Raw single-base T>C and A>G counts kept for QC. Not an rNMP signal "
            "and not shown as one."
        ),
    }

    _header, strand_rows = _table(strand_tsv)
    biases = []
    for row in strand_rows:
        if "strand_bias" in row and row["strand_bias"] != "":
            biases.append(float(row["strand_bias"]))
    avg_bias = round(sum(biases) / len(biases), 4) if biases else None

    with open(output_tsv, "w") as handle:
        handle.write("metric\tvalue\n")
        lines = [
            "headline_metric\trepeat_deletion_2_5bp\n",
            f"total_repeat_deletions\t{len(deletions)}\n",
        ]
        lines.extend(f"deletion_{length}bp\t{deletion_sizes.get(length, 0)}\n" for length in range(2, 6))
        lines.extend(f"repeat_unit_{unit}bp\t{unit_sizes.get(unit, 0)}\n" for unit in range(2, 6))
        handle.writelines(lines)

    signature = {
        "assay_type": "wgs",
        "headline_metric": "repeat_deletion_2_5bp",
        "headline_note": REPEAT_DELETION_NOTE,
        "total_repeat_deletions": len(deletions),
        "deletion_size_distribution": {str(length): deletion_sizes.get(length, 0) for length in range(2, 6)},
        "repeat_unit_distribution": {str(unit): unit_sizes.get(unit, 0) for unit in range(2, 6)},
        "avg_strand_bias": avg_bias,
        "polymerase_assignment": "experimental_not_assigned",
        "polymerase_note": POLYMERASE_NOTE,
        "v2_note": V2_POLYMERASE_NOTE,
        "accession": accession or "",
        "genotype": genotype or "",
        "qc": qc,
    }
    Path(signature_json).write_text(json.dumps(signature, indent=2) + "\n")
    print(
        f"WGS scan: {len(deletions)} repeat deletions (2–5 bp) → {output_tsv}",
        file=sys.stderr,
    )
    return signature


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize an rNMP signal table")
    parser.add_argument("--signal", required=True)
    parser.add_argument("--strand", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--assay", default="ribose_seq")
    parser.add_argument("--technique", default=None)
    parser.add_argument("--fai", default=None)
    parser.add_argument("--stats", default=None)
    parser.add_argument("--bam", default=None, help="Accepted for compatibility; not read")
    parser.add_argument("--accession", default=None)
    parser.add_argument("--genotype", default=None)
    args = parser.parse_args()
    if args.bam:
        print("--bam is ignored", file=sys.stderr)

    if args.assay == "wgs":
        scan_wgs(
            args.signal,
            args.strand,
            args.output,
            args.signature,
            stats_json=args.stats,
            accession=args.accession,
            genotype=args.genotype,
        )
        return
    genome_bp = genome_bp_from_fai(args.fai) if args.fai else None
    scan_direct(
        args.signal,
        args.strand,
        args.output,
        args.signature,
        args.assay,
        args.technique,
        genome_bp,
        args.stats,
        accession=args.accession,
        genotype=args.genotype,
    )


if __name__ == "__main__":
    main()
