#!/usr/bin/env python3
"""Summarize rNMP calls.

Window classes are descriptive count ratios. They are not polymerase
assignments. The WGS T>C / A>G tally is kept and marked experimental.
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


def genome_bp_from_fai(path: str) -> int:
    total = 0
    with open(path) as handle:
        for line in handle:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1].isdigit():
                total += int(parts[1])
    return total


def scan_direct(signal_tsv, asymmetry_tsv, output_tsv, signature_json, assay, technique, genome_bp, stats_json):
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
        "polymerase_note": (
            "Window classes are descriptive ratios of rNMP counts "
            "(forward / (forward + reverse)). They are not Pol epsilon versus "
            "Pol alpha/delta assignments; replication origin and fork direction are not modeled."
        ),
        "call_stats": stats,
    }
    Path(signature_json).write_text(json.dumps(signature, indent=2) + "\n")
    print(
        f"{technique or assay}: {total_ends} rNMP counts, "
        f"{len(signal_rows)} positions, bases {dict(base_counts)} → {output_tsv}",
        file=sys.stderr,
    )
    return signature


def scan_wgs(signal_tsv, strand_tsv, output_tsv, signature_json):
    transitions = []
    deletions = []
    with open(signal_tsv) as handle:
        handle.readline()
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, mtype = parts[:5]
            if mtype == "rNMP_transition":
                transitions.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt})
            elif mtype == "rNMP_deletion_slippage":
                del_len = int(parts[8]) if len(parts) > 8 and parts[8] else 0
                deletions.append({"chrom": chrom, "pos": int(pos), "del_len": del_len})

    transition_types = Counter(f"{row['ref']}→{row['alt']}" for row in transitions)
    deletion_sizes = Counter(row["del_len"] for row in deletions)
    total = len(transitions) + len(deletions)

    _header, strand_rows = _table(strand_tsv)
    biases = []
    for row in strand_rows:
        if "strand_bias" in row:
            biases.append(float(row["strand_bias"]))
    avg_bias = sum(biases) / len(biases) if biases else 0.0

    with open(output_tsv, "w") as handle:
        handle.write("type\tcount\tfrac_of_total\n")
        lines = [
            f"transition_{name}\t{count}\t{count / max(total, 1):.4f}\n"
            for name, count in transition_types.items()
        ]
        lines.extend(
            f"deletion_{length}bp\t{count}\t{count / max(total, 1):.4f}\n"
            for length, count in sorted(deletion_sizes.items())
        )
        lines.append(f"TOTAL\t{total}\t1.0\n")
        handle.writelines(lines)

    signature = {
        "assay_type": "wgs",
        "experimental": True,
        "interpretation": (
            "T>C and A>G counts are an experimental heuristic retained for review. "
            "They are not a validated rNMP mutation signature."
        ),
        "total_rnmp_transitions": len(transitions),
        "total_rnmp_deletions": len(deletions),
        "transition_breakdown": dict(transition_types),
        "deletion_size_distribution": {str(key): value for key, value in deletion_sizes.items()},
        "avg_strand_bias": round(avg_bias, 4),
        "rnmp_enrichment_score": round(len(transitions) / max(total, 1), 4),
        "polymerase_assignment": "experimental_not_assigned",
    }
    Path(signature_json).write_text(json.dumps(signature, indent=2) + "\n")
    print(
        f"WGS scan (experimental): {len(transitions)} transitions, {len(deletions)} deletions",
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
    args = parser.parse_args()
    if args.bam:
        print("--bam is ignored", file=sys.stderr)

    if args.assay == "wgs":
        scan_wgs(args.signal, args.strand, args.output, args.signature)
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
    )


if __name__ == "__main__":
    main()
