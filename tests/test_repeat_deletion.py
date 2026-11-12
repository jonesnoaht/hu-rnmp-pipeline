"""2–5 bp short-tandem-repeat deletions are the headline signature.

T>C and A>G stay under the JSON qc key and are not written as signal rows
or rendered in the report.
"""

import pysam
from build_report import build_report
from rnmp_mutation_scan import scan_wgs
from rnmp_signal import (
    classify_repeat_deletion,
    extract_wgs_rnmp_signal,
    simple_deletion,
    tandem_unit_covering,
)


def test_str_detection_accepts_di_and_tri_repeats():
    # One AT unit needs the neighboring copies to be a tandem repeat.
    assert classify_repeat_deletion("GAT", "G", left_flank="ATAT", right_flank="ATAT") == {
        "type": "repeat_deletion",
        "deleted": "AT",
        "del_len": 2,
        "repeat_unit": 2,
    }
    # The allele itself is already (AT)x2.
    assert classify_repeat_deletion("GATAT", "G")["repeat_unit"] == 2
    assert classify_repeat_deletion("GATG", "G", left_flank="ATG", right_flank="ATG")["repeat_unit"] == 3
    # 5 bp sits inside (AT)n even though 5 is not a multiple of 2.
    five = classify_repeat_deletion("GATATA", "G", left_flank="AT", right_flank="ATAT")
    assert five["del_len"] == 5
    assert five["repeat_unit"] == 2
    # One copy of a tetranucleotide unit, with a second copy in the flank.
    tetra = classify_repeat_deletion("GACGT", "G", left_flank="", right_flank="ACGT")
    assert tetra["del_len"] == 4
    assert tetra["repeat_unit"] == 4


def test_str_detection_rejects_unique_sequence_and_homopolymers():
    assert classify_repeat_deletion("GTCG", "G", left_flank="AA", right_flank="CC") is None
    assert classify_repeat_deletion("GAT", "G") is None  # a single AT is not a tandem repeat
    # Four identical bases are a homopolymer, not a di-nucleotide repeat.
    assert classify_repeat_deletion("GAAAA", "G", left_flank="AAAA", right_flank="AAAA") is None
    assert tandem_unit_covering("AAAAAAAA", 2, 6) is None
    # The deletion has to be covered by the repeat, not merely adjacent to it.
    assert classify_repeat_deletion("GGAT", "G", left_flank="AT", right_flank="ATAT") is None


def test_deletion_bounds_are_2_to_5_bp():
    flanks = {"left_flank": "ATATATAT", "right_flank": "ATATATAT"}
    assert classify_repeat_deletion("GA", "G", **flanks) is None
    assert simple_deletion("GA", "G")[1] == "A"
    six = "G" + ("AT" * 3)
    assert simple_deletion(six, "G")[1] == "ATATAT"
    assert classify_repeat_deletion(six, "G", **flanks) is None
    for length in (2, 3, 4, 5):
        deleted = ("AT" * 3)[:length]
        called = classify_repeat_deletion("G" + deleted, "G", **flanks)
        assert called is not None
        assert called["del_len"] == length


def _reference() -> str:
    # 0: T>C. 1-2: GG. 3-12: (AT)x5. 13-14: CC. 15-22: homopolymer A. 23: G for an A>G site at 22.
    return "T" + "GG" + "ATATATATAT" + "CC" + "AAAAAAAA" + "G"


def _write_fasta(path, sequence: str) -> None:
    path.write_text(f">chrI\n{sequence}\n")
    pysam.faidx(str(path))


def _write_vcf(path, rows: list[tuple[int, str, str]]) -> None:
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrI,length=80>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    body = "".join(f"chrI\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n" for pos, ref, alt in rows)
    path.write_text(header + body)


def _write_bam(path) -> None:
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrI", "LN": 80}]})
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        del bam


def _left_aligned(sequence: str, start0: int, length: int) -> tuple[int, str, str]:
    anchor = start0 - 1
    return anchor + 1, sequence[anchor : start0 + length], sequence[anchor]


def test_wgs_counter_writes_repeat_deletions_and_hides_transitions(tmp_path):
    sequence = _reference()
    assert sequence[0] == "T"
    assert sequence[22] == "A"
    fasta = tmp_path / "ref.fa"
    vcf = tmp_path / "calls.vcf"
    bam = tmp_path / "sample.bam"
    out = tmp_path / "rnmp_signal.tsv"
    _write_fasta(fasta, sequence)
    _write_bam(bam)
    two = _left_aligned(sequence, 5, 2)  # AT inside (AT)n
    five = _left_aligned(sequence, 5, 5)  # 5 bp inside (AT)n
    one = _left_aligned(sequence, 3, 1)  # 1 bp, below the bound
    six = _left_aligned(sequence, 3, 6)  # 6 bp, above the bound
    homo = _left_aligned(sequence, 16, 4)  # AAAA
    _write_vcf(
        vcf,
        [
            (1, "T", "C"),
            (23, "A", "G"),
            two,
            five,
            one,
            six,
            homo,
        ],
    )

    stats = extract_wgs_rnmp_signal(str(vcf), str(bam), str(out), fasta_path=str(fasta))
    text = out.read_text()
    assert stats["headline_metric"] == "repeat_deletion_2_5bp"
    assert stats["total_repeat_deletions"] == 2
    assert stats["qc"]["raw_tc_ag_counts"] == {"T>C": 1, "A>G": 1, "total": 2}
    assert "repeat_deletion" in text
    assert "\t2\t2\tAT\n" in text
    assert "\t5\t2\tATATA\n" in text
    assert "T\tC" not in text
    assert "A\tG" not in text

    # Without reference flanks, a single AT allele is not itself a tandem repeat.
    bare = extract_wgs_rnmp_signal(str(vcf), str(bam), str(tmp_path / "nofasta.tsv"))
    assert bare["total_repeat_deletions"] == 0
    assert bare["qc"]["raw_tc_ag_counts"]["total"] == 2


def test_report_does_not_show_transition_qc_as_signal(tmp_path):
    signal = tmp_path / "signal.tsv"
    strand = tmp_path / "strand.tsv"
    mutations = tmp_path / "mut.tsv"
    signature = tmp_path / "sig.json"
    signal.write_text(
        "chrom\tpos\tref\talt\ttype\tdel_len\trepeat_unit\tdeleted\n"
        "chrI\t20\tGATAT\tG\trepeat_deletion\t4\t2\tATAT\n"
    )
    strand.write_text("chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n")
    stats = tmp_path / "stats.json"
    stats.write_text(
        '{"qc": {"raw_tc_ag_counts": {"T>C": 3, "A>G": 1, "total": 4}}}\n'
    )
    payload = scan_wgs(str(signal), str(strand), str(mutations), str(signature), stats_json=str(stats))
    assert payload["qc"]["raw_tc_ag_counts"]["total"] == 4
    assert "T>C" not in mutations.read_text()

    for name in ("enrich.tsv", "hits.json", "harm.tsv", "harm.json", "qc"):
        path = tmp_path / name
        if name.endswith("json"):
            path.write_text("{}\n")
        elif name == "qc":
            path.mkdir()
        else:
            path.write_text("metric\tvalue\n")
    report = tmp_path / "report"
    build_report(
        "SRR11364929",
        str(tmp_path / "qc"),
        str(signal),
        str(mutations),
        str(signature),
        str(tmp_path / "enrich.tsv"),
        str(tmp_path / "hits.json"),
        str(tmp_path / "harm.tsv"),
        str(tmp_path / "harm.json"),
        str(report),
    )
    page = (report / "report.html").read_text()
    before, future = page.split('id="future-work"', 1)
    assert "Repeat deletions: <strong>1</strong>" in before
    assert "T>C" not in page
    assert "A>G" not in page
    assert "hydroxyurea" not in before.lower()
    assert "experimental_not_assigned" in before
    assert "Pol ε" not in page
    assert "Pol alpha" not in page
    assert "30847483" in future
    assert "25664722" in future
    assert "hydroxyurea" in future.lower()
    saved = (report / "mutation_signature.json").read_text()
    assert '"qc"' in saved
    assert '"T>C": 3' in saved
