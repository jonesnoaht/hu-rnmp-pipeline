#!/usr/bin/env python3
"""Extract rNMP positions from an alignment, with one rule per assay.

Direct assays use the Ribose-Map coordinate table in ``rnmp_coordinates``.
Only read 1 of a proper pair is counted. WGS mode still records T>C and A>G
calls as an experimental heuristic; that classification is unchanged and is
labeled experimental by the mutation scan.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam
from rnmp_coordinates import (
    AlignmentView,
    canonical_technique,
    five_prime_clip_len,
    keep_for_coordinate,
    rnmp_interval,
)

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def incorporated_base(reference_base: str, rnmp_strand: str) -> str:
    base = (reference_base or "N").upper()[:1]
    if base not in "ACGTN":
        base = "N"
    if rnmp_strand == "-":
        return base.translate(_COMPLEMENT)
    return base


def view_from_read(read: pysam.AlignedSegment) -> AlignmentView:
    start = read.reference_start if read.reference_start is not None else -1
    end = read.reference_end if read.reference_end is not None else -1
    return AlignmentView(
        strand="-" if read.is_reverse else "+",
        start=start,
        end=end if end is not None else -1,
        mapq=int(read.mapping_quality),
        is_paired=bool(read.is_paired),
        is_read1=bool(read.is_read1),
        is_proper_pair=bool(read.is_proper_pair),
        is_secondary=bool(read.is_secondary),
        is_supplementary=bool(read.is_supplementary),
        is_unmapped=bool(read.is_unmapped),
        is_qcfail=bool(read.is_qcfail),
        five_prime_clip=five_prime_clip_len(read.cigartuples, read.is_reverse),
    )


def _window_rows(counts: Counter, window: int) -> list[dict]:
    buckets: dict[tuple[str, int], dict[str, int]] = defaultdict(lambda: {"+": 0, "-": 0})
    for (chrom, pos, strand), count in counts.items():
        start = (pos // window) * window
        buckets[(chrom, start)][strand] += count
    rows = []
    for (chrom, start) in sorted(buckets):
        fwd = buckets[(chrom, start)]["+"]
        rev = buckets[(chrom, start)]["-"]
        total = fwd + rev
        ratio = fwd / total if total else 0.0
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": start + window,
                "fwd": fwd,
                "rev": rev,
                "total": total,
                "ratio": ratio,
            }
        )
    return rows


def classify_bias(bias: float) -> str:
    """Descriptive label only. Not a Pol epsilon / alpha / delta call."""
    if 0.45 <= bias <= 0.55:
        return "symmetric"
    if bias > 0.65:
        return "forward_biased"
    if bias < 0.35:
        return "reverse_biased"
    return "mild_bias"


def asymmetry_rows(counts: Counter, chrom_lengths: dict[str, int], window: int, min_total: int) -> list[dict]:
    buckets: dict[tuple[str, int], dict[str, int]] = defaultdict(lambda: {"+": 0, "-": 0})
    for (chrom, pos, strand), count in counts.items():
        start = (pos // window) * window
        buckets[(chrom, start)][strand] += count
    rows = []
    for (chrom, start) in sorted(buckets):
        fwd = buckets[(chrom, start)]["+"]
        rev = buckets[(chrom, start)]["-"]
        total = fwd + rev
        if total < min_total:
            continue
        end = start + window
        length = chrom_lengths.get(chrom)
        if length is not None:
            end = min(end, length)
        bias = fwd / total
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "fwd_rnmp": fwd,
                "rev_rnmp": rev,
                "strand_bias": round(bias, 4),
                "asymmetry": classify_bias(bias),
            }
        )
    return rows


def extract_direct_signal(
    bam_path: str,
    output_tsv: str,
    technique: str,
    fasta_path: str | None = None,
    min_mapq: int = 0,
    window: int = 100,
    asymmetry_window: int = 10_000,
) -> dict:
    tech = canonical_technique(technique)
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    skipped: Counter[str] = Counter()
    counts: Counter[tuple[str, int, str]] = Counter()
    bases: dict[tuple[str, int, str], str] = {}
    kept = 0
    seen = 0
    off_chromosome = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_lengths = dict(zip(bam.references, bam.lengths))
        for read in bam.fetch(until_eof=True):
            seen += 1
            view = view_from_read(read)
            reason = keep_for_coordinate(view, min_mapq)
            if reason:
                skipped[reason] += 1
                continue
            chrom = bam.get_reference_name(read.reference_id)
            interval = rnmp_interval(
                tech,
                view.start,
                view.end,
                view.strand,
                chrom_lengths.get(chrom),
            )
            if interval is None:
                off_chromosome += 1
                continue
            start, _end, strand = interval
            key = (chrom, start, strand)
            counts[key] += 1
            kept += 1
            if fasta is not None and key not in bases:
                try:
                    ref_base = fasta.fetch(chrom, start, start + 1)
                except (ValueError, KeyError):
                    ref_base = "N"
                bases[key] = incorporated_base(ref_base, strand)

    if fasta is not None:
        fasta.close()

    with open(output_tsv, "w") as handle:
        handle.write("chrom\tpos\tstrand\tcount\tbase\n")
        handle.writelines(
            f"{chrom}\t{pos}\t{strand}\t{counts[(chrom, pos, strand)]}\t{bases.get((chrom, pos, strand), 'NA')}\n"
            for chrom, pos, strand in sorted(counts)
        )

    windows_path = output_tsv.replace(".tsv", "_windows.tsv")
    with open(windows_path, "w") as handle:
        handle.write("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tstrand_ratio\n")
        handle.writelines(
            f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['fwd']}\t{row['rev']}\t"
            f"{row['total']}\t{row['ratio']:.4f}\n"
            for row in _window_rows(counts, window)
        )

    asymmetry_path = output_tsv.replace(".tsv", "_asymmetry.tsv")
    with open(asymmetry_path, "w") as handle:
        handle.write("chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n")
        handle.writelines(
            f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['fwd_rnmp']}\t{row['rev_rnmp']}\t"
            f"{row['strand_bias']}\t{row['asymmetry']}\n"
            for row in asymmetry_rows(counts, chrom_lengths, asymmetry_window, min_total=5)
        )

    stats = {
        "technique": tech,
        "min_mapq": min_mapq,
        "reads_seen": seen,
        "reads_kept": kept,
        "off_chromosome": off_chromosome,
        "skipped": dict(skipped),
        "rnmp_sites": len(counts),
        "rnmp_counts": int(sum(counts.values())),
        "coordinate_reference": "Ribose-Map modules/coordinate.sh @ fff581bf",
        "windows": windows_path,
        "asymmetry": asymmetry_path,
    }
    stats_path = output_tsv.replace(".tsv", "_stats.json")
    Path(stats_path).write_text(json.dumps(stats, indent=2) + "\n")
    print(
        f"{tech}: kept {kept} of {seen} reads, {len(counts)} rNMP sites → {output_tsv}",
        file=sys.stderr,
    )
    return stats


def extract_wgs_rnmp_signal(vcf_path: str, bam_path: str, output_tsv: str) -> dict:
    """Experimental T>C / A>G counter. The deletion branch below is unreachable.

    ``len(ref) > 1`` continues before the 2-5 bp deletion test. That behavior
    is intentionally unchanged pending a decision on whether this signature stays.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    rnmp_transitions = {("T", "C"), ("A", "G")}
    signals = []

    try:
        vcf = pysam.VariantFile(vcf_path)
    except (OSError, ValueError):
        vcf = None

    if vcf:
        for rec in vcf.fetch():
            ref = rec.ref.upper()
            for alt in rec.alts or []:
                alt = alt.upper()
                if alt == "<NON_REF>" or len(alt) > 1 or len(ref) > 1:
                    continue
                if (ref, alt) in rnmp_transitions:
                    depth_fwd = 0
                    depth_rev = 0
                    for pileup in bam.pileup(rec.chrom, rec.pos - 1, rec.pos):
                        if pileup.pos == rec.pos - 1:
                            for read in pileup.pileups:
                                if read.is_forward:
                                    depth_fwd += 1
                                else:
                                    depth_rev += 1
                    signals.append(
                        {
                            "chrom": rec.chrom,
                            "pos": rec.pos,
                            "ref": ref,
                            "alt": alt,
                            "type": "rNMP_transition",
                            "depth_fwd": depth_fwd,
                            "depth_rev": depth_rev,
                            "strand_bias": depth_fwd / max(depth_fwd + depth_rev, 1),
                        }
                    )
                elif 2 <= len(ref) <= 5 and len(alt) == 1:
                    signals.append(
                        {
                            "chrom": rec.chrom,
                            "pos": rec.pos,
                            "ref": ref,
                            "alt": alt,
                            "type": "rNMP_deletion_slippage",
                            "del_len": len(ref) - 1,
                        }
                    )
    bam.close()

    columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "type",
        "depth_fwd",
        "depth_rev",
        "strand_bias",
        "del_len",
    ]
    with open(output_tsv, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.writelines("\t".join(str(row.get(key, "")) for key in columns) + "\n" for row in signals)

    print(f"WGS rNMP signal (experimental): {len(signals)} candidate sites → {output_tsv}", file=sys.stderr)
    return {"total_signals": len(signals), "experimental": True}


def write_empty_companions(output_tsv: str) -> None:
    windows = output_tsv.replace(".tsv", "_windows.tsv")
    asymmetry = output_tsv.replace(".tsv", "_asymmetry.tsv")
    stats = output_tsv.replace(".tsv", "_stats.json")
    Path(windows).write_text("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tstrand_ratio\n")
    Path(asymmetry).write_text("chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n")
    if not Path(stats).exists():
        Path(stats).write_text(json.dumps({"assay": "wgs", "experimental": True}) + "\n")


def identify_rnmp_hotspots(signal_tsv: str, output_tsv: str, min_count: int = 5, window: int = 500) -> list[dict]:
    """Windows at least 2x the median non-empty window. This threshold is arbitrary."""
    counts: dict[str, dict[tuple[int, str], int]] = defaultdict(lambda: defaultdict(int))
    with open(signal_tsv) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: pos for pos, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if "strand" not in index or len(parts) <= index.get("count", 3):
                continue
            try:
                chrom = parts[index.get("chrom", 0)]
                pos = int(parts[index.get("pos", 1)])
                strand = parts[index["strand"]]
                count = int(parts[index["count"]])
            except (ValueError, IndexError):
                continue
            counts[chrom][(pos, strand)] += count

    window_counts: dict[tuple[str, int], dict[str, int]] = defaultdict(lambda: {"+": 0, "-": 0})
    for chrom, pos_strands in counts.items():
        for (pos, strand), count in pos_strands.items():
            win_start = (pos // window) * window
            if strand in {"+", "-"}:
                window_counts[(chrom, win_start)][strand] += count

    totals = [strands["+"] + strands["-"] for strands in window_counts.values()]
    with open(output_tsv, "w") as handle:
        handle.write("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tfold_over_median\n")
        if not totals:
            print("No rNMP sites found for hotspot analysis", file=sys.stderr)
            return []
        totals.sort()
        median = totals[len(totals) // 2]
        hotspots = []
        for (chrom, start), strands in window_counts.items():
            fwd = strands["+"]
            rev = strands["-"]
            total = fwd + rev
            fold = total / median if median > 0 else 0
            if total >= min_count and fold >= 2.0:
                hotspots.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": start + window,
                        "fwd_count": fwd,
                        "rev_count": rev,
                        "total": total,
                        "fold_over_median": round(fold, 2),
                    }
                )
        handle.writelines(
            f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['fwd_count']}\t{row['rev_count']}\t"
            f"{row['total']}\t{row['fold_over_median']}\n"
            for row in sorted(hotspots, key=lambda row: row["total"], reverse=True)
        )
    print(
        f"Hotspots: {len(hotspots)} regions ≥2x median (median={median}) → {output_tsv}",
        file=sys.stderr,
    )
    return hotspots


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract rNMP signal from an alignment or VCF")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--assay", default="ribose_seq")
    parser.add_argument("--technique", default=None, help="Ribose-Map technique name")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--vcf", default=None)
    parser.add_argument("--hotspots", default=None)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--window", type=int, default=100)
    args = parser.parse_args()

    if args.assay == "wgs":
        if not args.vcf:
            parser.error("--vcf is required for --assay wgs")
        extract_wgs_rnmp_signal(args.vcf, args.bam, args.output)
        write_empty_companions(args.output)
    else:
        if not args.technique:
            parser.error("--technique is required for direct rNMP assays")
        extract_direct_signal(
            args.bam,
            args.output,
            args.technique,
            fasta_path=args.fasta,
            min_mapq=args.min_mapq,
            window=args.window,
        )
        if args.hotspots:
            identify_rnmp_hotspots(args.output, args.hotspots)


if __name__ == "__main__":
    main()
