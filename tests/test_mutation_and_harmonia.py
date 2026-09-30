"""Summary labels stay experimental, and Harmonia skips without a binary."""

import json
import sys

from harmonia_join import main as harmonia_main
from rnmp_mutation_scan import scan_direct, scan_wgs


def test_direct_scan_does_not_assign_polymerase(tmp_path):
    signal = tmp_path / "signal.tsv"
    asym = tmp_path / "asym.tsv"
    out = tmp_path / "mut.tsv"
    signature = tmp_path / "sig.json"
    signal.write_text(
        "chrom\tpos\tstrand\tcount\tbase\n"
        "chrI\t10\t+\t4\tA\n"
        "chrI\t11\t-\t1\tC\n"
    )
    asym.write_text(
        "chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n"
        "chrI\t0\t10000\t4\t1\t0.8\tforward_biased\n"
    )
    payload = scan_direct(str(signal), str(asym), str(out), str(signature), "ribose_seq", "ribose-seq", 1_000_000, None)
    assert payload["polymerase_assignment"] == "experimental_not_assigned"
    assert payload["technique"] == "ribose-seq"
    assert payload["base_counts"] == {"A": 4, "C": 1}
    assert payload["forward_ends"] == 4
    assert payload["reverse_ends"] == 1
    assert payload["rnmp_calls_per_mb"] == 5.0
    assert payload["rnmp_calls_per_mb_note"]
    text = out.read_text()
    assert "\nchrom\t" not in text


def test_wgs_signature_is_marked_experimental(tmp_path):
    signal = tmp_path / "signal.tsv"
    strand = tmp_path / "strand.tsv"
    out = tmp_path / "mut.tsv"
    signature = tmp_path / "sig.json"
    signal.write_text("chrom\tpos\tref\talt\ttype\nchrI\t10\tT\tC\trNMP_transition\n")
    strand.write_text("chrom\tstart\tend\tfwd_rnmp\trev_rnmp\tstrand_bias\tasymmetry\n")
    payload = scan_wgs(str(signal), str(strand), str(out), str(signature))
    assert payload["experimental"] is True
    assert payload["total_rnmp_transitions"] == 1
    assert payload["polymerase_assignment"] == "experimental_not_assigned"


def test_harmonia_skips_when_disabled(tmp_path, capsys):
    enrichment = tmp_path / "enrich.tsv"
    rules = tmp_path / "rules.txt"
    output = tmp_path / "matrix.tsv"
    complement = tmp_path / "complement.json"
    enrichment.write_text("status\tok\n")
    rules.write_text("rule\n")
    argv = sys.argv
    sys.argv = [
        "harmonia_join.py",
        "--enrichment",
        str(enrichment),
        "--rules",
        str(rules),
        "--output",
        str(output),
        "--complement",
        str(complement),
        "--run",
        "false",
    ]
    try:
        harmonia_main()
    finally:
        sys.argv = argv
    assert "skipped" in output.read_text()
    assert json.loads(complement.read_text())["status"] == "skipped"
    assert capsys.readouterr().err
