"""Coordinate rules from Ribose-Map coordinate.sh @ fff581b (PMID 34089018 Table 1)."""

import pytest
from rnmp_coordinates import (
    AlignmentView,
    canonical_technique,
    keep_for_coordinate,
    rnmp_interval,
)

READ_START = 100
READ_END = 150


@pytest.mark.parametrize(
    ("technique", "strand", "expected"),
    [
        ("ribose-seq", "+", (100, 101, "-")),
        ("ribose-seq", "-", (149, 150, "+")),
        ("emRiboSeq", "+", (99, 100, "-")),
        ("emRiboSeq", "-", (150, 151, "+")),
        ("Alk-HydEn-seq", "+", (99, 100, "+")),
        ("Alk-HydEn-seq", "-", (150, 151, "-")),
        ("Pu-seq", "+", (99, 100, "+")),
        ("Pu-seq", "-", (150, 151, "-")),
        ("RHII-HydEn-seq", "+", (100, 101, "+")),
        ("RHII-HydEn-seq", "-", (149, 150, "-")),
    ],
)
def test_coordinate_table(technique, strand, expected):
    assert rnmp_interval(technique, READ_START, READ_END, strand, chrom_length=1000) == expected


def test_aliases_match_canonical_rules():
    assert canonical_technique("ribose_seq") == "ribose-seq"
    assert rnmp_interval("alk_hyden_seq", 100, 150, "+", 1000) == rnmp_interval(
        "Alk-HydEn-seq", 100, 150, "+", 1000
    )
    assert rnmp_interval("rhii_hyden_seq", 100, 150, "-", 1000) == (149, 150, "-")


def test_bare_hyden_seq_is_rejected():
    with pytest.raises(ValueError, match="ambiguous"):
        rnmp_interval("hyden_seq", 100, 150, "+", 1000)


def test_alkaline_upstream_of_chromosome_start_is_dropped():
    assert rnmp_interval("Alk-HydEn-seq", 0, 50, "+", chrom_length=1000) is None
    assert rnmp_interval("emRiboSeq", 0, 50, "+", chrom_length=1000) is None


def test_minus_strand_past_chromosome_end_is_dropped():
    assert rnmp_interval("Alk-HydEn-seq", 950, 1000, "-", chrom_length=1000) is None
    assert rnmp_interval("Pu-seq", 950, 1000, "-", chrom_length=1000) is None
    assert rnmp_interval("emRiboSeq", 950, 1000, "-", chrom_length=1000) is None


def _view(**overrides) -> AlignmentView:
    fields = {
        "strand": "+",
        "start": 100,
        "end": 150,
        "mapq": 40,
        "is_paired": True,
        "is_read1": True,
        "is_proper_pair": True,
        "is_secondary": False,
        "is_supplementary": False,
        "is_unmapped": False,
        "is_qcfail": False,
        "five_prime_clip": 0,
    }
    fields.update(overrides)
    return AlignmentView(**fields)


def test_paired_end_keeps_read1_only():
    assert keep_for_coordinate(_view(), min_mapq=0) is None
    assert keep_for_coordinate(_view(is_read1=False), min_mapq=0) == "not_read1_proper_pair"
    assert keep_for_coordinate(_view(is_proper_pair=False), min_mapq=0) == "not_read1_proper_pair"


def test_single_end_is_kept():
    assert keep_for_coordinate(_view(is_paired=False, is_read1=False, is_proper_pair=False), min_mapq=0) is None


def test_mapq_and_five_prime_clip_filters():
    assert keep_for_coordinate(_view(mapq=0), min_mapq=1) == "mapq"
    assert keep_for_coordinate(_view(mapq=0), min_mapq=0) is None
    assert keep_for_coordinate(_view(five_prime_clip=3), min_mapq=0) == "five_prime_clip"
