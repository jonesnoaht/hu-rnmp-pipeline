#!/usr/bin/env python3
"""Extract rNMP positions from an alignment, with one rule per assay.

Direct assays use the Ribose-Map coordinate table in ``rnmp_coordinates``.
Only read 1 of a proper pair is counted. WGS mode's headline mutation
signature is a 2–5 bp deletion that sits in a short tandem repeat. Raw
single-base T>C and A>G counts are returned under ``qc`` and are not an
rNMP signal.
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


# Deletion must be this long, and the repeat unit is a di-/tri-/tetra-/pentanucleotide.
# Homopolymer units are excluded: a run of identical bases is not the short-tandem-repeat
# signature described for Top1 cleavage at an embedded rNMP.
REPEAT_DELETION_MIN = 2
REPEAT_DELETION_MAX = 5
REPEAT_UNIT_LENGTHS = (2, 3, 4, 5)
REPEAT_FLANK = 12
_TRANSITION_QC = {("T", "C"): "T>C", ("A", "G"): "A>G"}


def classify_bias(bias: float) -> str:
    """Descriptive strand-count label. v1 does not assign replicative polymerases."""
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


def _is_homopolymer(unit: str) -> bool:
    return bool(unit) and all(base == unit[0] for base in unit)


def simple_deletion(ref: str, alt: str, pos: int = 1) -> tuple[int, str] | None:
    """Return the 0-based genomic start and deleted sequence of a simple deletion.

    Left-aligned alleles (``REF`` starts with ``ALT``) and right-aligned alleles
    (``REF`` ends with ``ALT``) are both accepted. Anything that is not a pure
    deletion of ACGT bases returns ``None``. ``pos`` is the 1-based VCF POS.
    """
    ref_u = (ref or "").upper()
    alt_u = (alt or "").upper()
    if not ref_u or not alt_u:
        return None
    if any(base not in "ACGT" for base in ref_u + alt_u):
        return None
    if len(ref_u) <= len(alt_u):
        return None
    if ref_u.startswith(alt_u):
        deleted = ref_u[len(alt_u) :]
        start0 = (pos - 1) + len(alt_u)
        return start0, deleted
    if ref_u.endswith(alt_u):
        deleted = ref_u[: len(ref_u) - len(alt_u)]
        return pos - 1, deleted
    return None


def transition_qc_key(ref: str, alt: str) -> str | None:
    """Return ``T>C`` or ``A>G`` for a single-base change, else ``None``.

    These counts are a QC tally. They are not the rNMP mutation signature.
    """
    ref_u = (ref or "").upper()
    alt_u = (alt or "").upper()
    if len(ref_u) != 1 or len(alt_u) != 1:
        return None
    return _TRANSITION_QC.get((ref_u, alt_u))


def tandem_unit_covering(
    sequence: str,
    start: int,
    end: int,
    unit_lengths: tuple[int, ...] = REPEAT_UNIT_LENGTHS,
) -> int | None:
    """Unit length of a short tandem repeat that fully covers ``sequence[start:end]``.

    The repeat must be at least two copies of a 2–5 bp unit. A unit whose bases
    are all the same is a homopolymer and does not count. The shortest covering
    unit is returned.
    """
    sequence = (sequence or "").upper()
    if start < 0 or end > len(sequence) or start >= end:
        return None
    for unit_len in unit_lengths:
        for phase in range(unit_len):
            index = phase
            run_start = None
            previous = None
            while index + unit_len <= len(sequence):
                unit = sequence[index : index + unit_len]
                continues = previous is not None and unit == previous and not _is_homopolymer(unit)
                if continues:
                    index += unit_len
                    continue
                if (
                    run_start is not None
                    and previous is not None
                    and not _is_homopolymer(previous)
                    and run_start <= start
                    and index >= end
                    and (index - run_start) >= 2 * unit_len
                ):
                    return unit_len
                if _is_homopolymer(unit):
                    run_start = None
                    previous = None
                else:
                    run_start = index
                    previous = unit
                index += unit_len
            if (
                run_start is not None
                and previous is not None
                and not _is_homopolymer(previous)
                and run_start <= start
                and index >= end
                and (index - run_start) >= 2 * unit_len
            ):
                return unit_len
    return None


def classify_repeat_deletion(ref: str, alt: str, left_flank: str = "", right_flank: str = "") -> dict | None:
    """Classify a 2–5 bp deletion that sits inside a short tandem repeat.

    ``left_flank`` and ``right_flank`` are reference bases immediately before
    and after the deleted sequence. A deletion whose own allele is already two
    or more repeat units can qualify with empty flanks. Returns ``None`` when
    the length is outside 2–5 bp or the deletion is not inside a tandem repeat.
    """
    parsed = simple_deletion(ref, alt)
    if parsed is None:
        return None
    _start, deleted = parsed
    length = len(deleted)
    if length < REPEAT_DELETION_MIN or length > REPEAT_DELETION_MAX:
        return None
    window = f"{left_flank}{deleted}{right_flank}".upper()
    unit = tandem_unit_covering(window, len(left_flank), len(left_flank) + length)
    if unit is None:
        return None
    return {
        "type": "repeat_deletion",
        "deleted": deleted,
        "del_len": length,
        "repeat_unit": unit,
    }


def _reference_flanks(fasta, chrom: str, start0: int, length: int) -> tuple[str, str]:
    if fasta is None:
        return "", ""
    left_start = max(0, start0 - REPEAT_FLANK)
    try:
        left = fasta.fetch(chrom, left_start, start0)
        right = fasta.fetch(chrom, start0 + length, start0 + length + REPEAT_FLANK)
    except (ValueError, KeyError, IndexError):
        return "", ""
    return left.upper(), right.upper()


def extract_wgs_rnmp_signal(
    vcf_path: str,
    bam_path: str,
    output_tsv: str,
    fasta_path: str | None = None,
) -> dict:
    """Count 2–5 bp repeat deletions. Keep T>C / A>G only under ``qc``.

    The deletion has to sit in a short tandem repeat of a 2–5 bp unit (at least
    two copies). Homopolymer runs are not counted. Single-base T>C and A>G
    calls are tallied for QC and are not written as rNMP signal rows.
    """
    # Opened so a missing alignment fails the stage. Depth is not part of the
    # repeat-deletion call.
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    deletions: list[dict] = []
    qc_counts: Counter[str] = Counter()
    seen = 0
    rejected_length = 0
    rejected_not_repeat = 0

    try:
        vcf = pysam.VariantFile(vcf_path)
    except (OSError, ValueError):
        vcf = None

    if vcf is not None:
        for rec in vcf:
            ref = (rec.ref or "").upper()
            for alt in rec.alts or []:
                alt_u = alt.upper()
                if alt_u in {"<NON_REF>", ".", "*"}:
                    continue
                seen += 1
                qc_key = transition_qc_key(ref, alt_u)
                if qc_key:
                    qc_counts[qc_key] += 1
                    continue
                parsed = simple_deletion(ref, alt_u, rec.pos)
                if parsed is None:
                    continue
                start0, deleted = parsed
                length = len(deleted)
                if length < REPEAT_DELETION_MIN or length > REPEAT_DELETION_MAX:
                    rejected_length += 1
                    continue
                left, right = _reference_flanks(fasta, rec.chrom, start0, length)
                classified = classify_repeat_deletion(ref, alt_u, left, right)
                if classified is None:
                    rejected_not_repeat += 1
                    continue
                deletions.append(
                    {
                        "chrom": rec.chrom,
                        "pos": start0 + 1,
                        "ref": ref,
                        "alt": alt_u,
                        "type": "repeat_deletion",
                        "del_len": classified["del_len"],
                        "repeat_unit": classified["repeat_unit"],
                        "deleted": classified["deleted"],
                    }
                )
        vcf.close()
    bam.close()
    if fasta is not None:
        fasta.close()

    columns = ["chrom", "pos", "ref", "alt", "type", "del_len", "repeat_unit", "deleted"]
    with open(output_tsv, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.writelines("\t".join(str(row[key]) for key in columns) + "\n" for row in deletions)

    qc = {
        "raw_tc_ag_counts": {
            "T>C": int(qc_counts["T>C"]),
            "A>G": int(qc_counts["A>G"]),
            "total": int(qc_counts["T>C"] + qc_counts["A>G"]),
        },
        "note": (
            "Raw single-base T>C and A>G counts kept for QC. Not an rNMP signal. "
            "Embedded rNMPs in yeast are associated with Top1-dependent 2–5 bp "
            "deletions in short tandem repeats (Nick McElhinny et al. 2010, "
            "PMID 20729855; Kim et al. 2011, PMID 21700875)."
        ),
    }
    stats = {
        "assay": "wgs",
        "headline_metric": "repeat_deletion_2_5bp",
        "total_repeat_deletions": len(deletions),
        "variants_seen": seen,
        "rejected_deletion_length": rejected_length,
        "rejected_not_short_tandem_repeat": rejected_not_repeat,
        "qc": qc,
    }
    stats_path = output_tsv.replace(".tsv", "_stats.json")
    Path(stats_path).write_text(json.dumps(stats, indent=2) + "\n")
    print(
        f"WGS repeat deletions: {len(deletions)} sites of 2–5 bp in short tandem repeats → {output_tsv}",
        file=sys.stderr,
    )
    return stats


def write_empty_companions(output_tsv: str) -> None:
    windows = output_tsv.replace(".tsv", "_windows.tsv")
    asymmetry = output_tsv.replace(".tsv", "_asymmetry.tsv")
    stats = output_tsv.replace(".tsv", "_stats.json")
    Path(windows).write_text("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tstrand_ratio\n")
    Path(asymmetry).write_text("chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n")
    if not Path(stats).exists():
        Path(stats).write_text(
            json.dumps(
                {
                    "assay": "wgs",
                    "headline_metric": "repeat_deletion_2_5bp",
                    "total_repeat_deletions": 0,
                    "qc": {"raw_tc_ag_counts": {"T>C": 0, "A>G": 0, "total": 0}},
                }
            )
            + "\n"
        )


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
        extract_wgs_rnmp_signal(args.vcf, args.bam, args.output, fasta_path=args.fasta)
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
