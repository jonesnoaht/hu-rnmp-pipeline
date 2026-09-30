"""BAM-level check that read 2 and 5' clips are not counted, and the assay rule is applied."""

import pysam
from rnmp_signal import extract_direct_signal


def _write_bam(path, reads):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrI", "LN": 1000}]}
    )
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for spec in reads:
            length = spec.get("length", 50)
            aln = pysam.AlignedSegment(header)
            aln.query_name = spec["name"]
            aln.query_sequence = "A" * length
            aln.query_qualities = pysam.qualitystring_to_array("I" * length)
            aln.reference_id = 0
            aln.reference_start = spec["start"]
            aln.cigar = spec.get("cigar", [(0, length)])
            aln.mapping_quality = spec.get("mapq", 40)
            aln.flag = spec["flag"]
            bam.write(aln)


def _positions(tsv):
    rows = []
    lines = tsv.read_text().splitlines()[1:]
    for line in lines:
        chrom, pos, strand, count, base = line.split("\t")
        rows.append((chrom, int(pos), strand, int(count), base))
    return rows


def test_rhii_counts_read1_not_read2(tmp_path):
    bam = tmp_path / "sample.bam"
    _write_bam(
        bam,
        [
            {"name": "r1", "start": 100, "flag": 67},
            {"name": "r2", "start": 400, "flag": 131},
            {"name": "r1rev", "start": 300, "flag": 83},
        ],
    )
    out = tmp_path / "rnmp_signal.tsv"
    extract_direct_signal(str(bam), str(out), "RHII-HydEn-seq")
    assert _positions(out) == [
        ("chrI", 100, "+", 1, "NA"),
        ("chrI", 349, "-", 1, "NA"),
    ]


def test_ribose_seq_flips_strand_and_alkaline_shifts(tmp_path):
    bam = tmp_path / "sample.bam"
    _write_bam(bam, [{"name": "se", "start": 100, "flag": 0}])
    ribose = tmp_path / "ribose.tsv"
    alk = tmp_path / "alk.tsv"
    extract_direct_signal(str(bam), str(ribose), "ribose-seq")
    extract_direct_signal(str(bam), str(alk), "Alk-HydEn-seq")
    assert _positions(ribose) == [("chrI", 100, "-", 1, "NA")]
    assert _positions(alk) == [("chrI", 99, "+", 1, "NA")]


def test_soft_clip_and_low_mapq_are_skipped(tmp_path):
    bam = tmp_path / "sample.bam"
    _write_bam(
        bam,
        [
            {"name": "clip", "start": 100, "flag": 0, "cigar": [(4, 5), (0, 45)]},
            {"name": "low", "start": 200, "flag": 0, "mapq": 0},
        ],
    )
    out = tmp_path / "rnmp_signal.tsv"
    stats = extract_direct_signal(str(bam), str(out), "ribose-seq", min_mapq=10)
    assert _positions(out) == []
    assert stats["skipped"]["five_prime_clip"] == 1
    assert stats["skipped"]["mapq"] == 1
