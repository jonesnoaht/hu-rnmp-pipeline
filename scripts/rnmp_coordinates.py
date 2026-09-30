#!/usr/bin/env python3
"""Assay-specific rNMP genomic coordinates.

The arithmetic matches Ribose-Map ``modules/coordinate.sh`` at commit
fff581bf28ff3e6c3970da99ccfb29df8c8de439 and Table 1 of Gombolay & Storici,
Nat Protoc 2021 (PMID 34089018). Intervals are 0-based, half-open, as in
``bedtools bamtobed``: read start S, read end E.

Paired-end data keeps read 1 of a proper pair (samtools ``-f 67``).
Single-end reads are kept. Alignments with a 5' soft or hard clip are
skipped, because those clips move the tagged end off ``reference_start``.
"""

from __future__ import annotations

from dataclasses import dataclass

# Canonical Ribose-Map technique names.
RIBOSE_SEQ = "ribose-seq"
EMRIBOSEQ = "emRiboSeq"
ALK_HYDEN = "Alk-HydEn-seq"
RHII_HYDEN = "RHII-HydEn-seq"
PU_SEQ = "Pu-seq"

_ALIASES = {
    "ribose-seq": RIBOSE_SEQ,
    "ribose_seq": RIBOSE_SEQ,
    "emriboseq": EMRIBOSEQ,
    "emRiboSeq": EMRIBOSEQ,
    "alk-hyden-seq": ALK_HYDEN,
    "alk_hyden_seq": ALK_HYDEN,
    "Alk-HydEn-seq": ALK_HYDEN,
    "pu-seq": PU_SEQ,
    "pu_seq": PU_SEQ,
    "Pu-seq": PU_SEQ,
    "rhii-hyden-seq": RHII_HYDEN,
    "rhii_hyden_seq": RHII_HYDEN,
    "RHII-HydEn-seq": RHII_HYDEN,
}

# Bare "hyden_seq" is not mapped: alkaline and RNase HII chemistries differ by 1 nt.
_AMBIGUOUS = {"hyden_seq", "hyden-seq", "hyden"}


def canonical_technique(name: str) -> str:
    """Return the Ribose-Map technique name, or raise ValueError."""
    key = (name or "").strip()
    if key in _AMBIGUOUS:
        raise ValueError(
            f"technique {key!r} is ambiguous. Use {RHII_HYDEN!r} "
            "(RNase HII nick; rNMP is the read 5' base, same strand) or "
            f"{ALK_HYDEN!r} (alkaline hydrolysis; rNMP is 1 nt upstream, same strand)."
        )
    try:
        return _ALIASES[key]
    except KeyError as exc:
        known = ", ".join(sorted(set(_ALIASES.values())))
        raise ValueError(f"Unknown technique {name!r}. Expected one of: {known}.") from exc


def rnmp_interval(
    technique: str,
    read_start: int,
    read_end: int,
    read_strand: str,
    chrom_length: int | None = None,
) -> tuple[int, int, str] | None:
    """Return ``(start, end, strand)`` of the rNMP, or None if it falls off the chromosome.

    ``read_start`` / ``read_end`` are the aligned read interval (0-based, half-open).
    ``read_strand`` is ``+`` or ``-`` as reported by ``bedtools bamtobed``.
    """
    tech = canonical_technique(technique)
    if read_strand not in {"+", "-"}:
        raise ValueError(f"read_strand must be '+' or '-', got {read_strand!r}")
    if read_end < read_start:
        raise ValueError(f"read interval is empty or reversed: {read_start}-{read_end}")

    if tech == RIBOSE_SEQ:
        if read_strand == "+":
            start, end, strand = read_start, read_start + 1, "-"
        else:
            start, end, strand = read_end - 1, read_end, "+"
    elif tech == EMRIBOSEQ:
        if read_strand == "+":
            start, end, strand = read_start - 1, read_start, "-"
        else:
            start, end, strand = read_end, read_end + 1, "+"
    elif tech in {ALK_HYDEN, PU_SEQ}:
        if read_strand == "+":
            start, end, strand = read_start - 1, read_start, "+"
        else:
            start, end, strand = read_end, read_end + 1, "-"
    elif tech == RHII_HYDEN:
        if read_strand == "+":
            start, end, strand = read_start, read_start + 1, "+"
        else:
            start, end, strand = read_end - 1, read_end, "-"
    else:
        raise ValueError(f"no coordinate rule for {tech}")

    if start < 0:
        return None
    if chrom_length is not None and end > chrom_length:
        return None
    return (start, end, strand)


@dataclass(frozen=True)
class AlignmentView:
    """The subset of a SAM record the coordinate filter needs."""

    strand: str
    start: int
    end: int
    mapq: int
    is_paired: bool
    is_read1: bool
    is_proper_pair: bool
    is_secondary: bool
    is_supplementary: bool
    is_unmapped: bool
    is_qcfail: bool
    five_prime_clip: int


def keep_for_coordinate(aln: AlignmentView, min_mapq: int) -> str | None:
    """Return None when the read carries an rNMP site, else a skip reason.

    Paired-end matches Ribose-Map ``samtools view -f 67`` (paired, proper pair,
    read 1). Secondary and supplementary alignments are also dropped so one
    fragment is not counted twice. Reads with a 5' clip are dropped rather
    than shifted: the tagged base is not at the aligned start.
    """
    if aln.is_unmapped:
        return "unmapped"
    if aln.is_secondary:
        return "secondary"
    if aln.is_supplementary:
        return "supplementary"
    if aln.is_qcfail:
        return "qcfail"
    if aln.mapq < min_mapq:
        return "mapq"
    if aln.five_prime_clip > 0:
        return "five_prime_clip"
    if aln.is_paired and not (aln.is_read1 and aln.is_proper_pair):
        return "not_read1_proper_pair"
    return None


def five_prime_clip_len(cigartuples: list[tuple[int, int]] | None, is_reverse: bool) -> int:
    """Length of a 5' soft clip (op 4) or hard clip (op 5), in query bases."""
    if not cigartuples:
        return 0
    op, length = cigartuples[-1] if is_reverse else cigartuples[0]
    if op in {4, 5}:
        return length
    return 0
