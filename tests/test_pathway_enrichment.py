"""The constant panel-size score is gone. Tests use the signal counts."""

from pathway_enrichment import (
    hypergeometric_pvalue,
    load_metric_table,
    run,
)


def _fasta(path):
    path.write_text(">chrI\n" + ("ACGT" * 50) + "\n")


def _signal(path, rows):
    lines = ["chrom\tpos\tstrand\tcount\tbase\n"]
    lines.extend(f"{chrom}\t{pos}\t{strand}\t{count}\t{base}\n" for chrom, pos, strand, count, base in rows)
    path.write_text("".join(lines))


def _pathways(path):
    path.write_text(
        "pathways:\n"
        "  - id: rnr-regulation\n"
        "    name: RNR\n"
        "    rationale: placeholder text\n"
    )


def _mutations(path, value="0.6250"):
    path.write_text(f"metric\tvalue\nforward_fraction\t{value}\n\nchrom\tfwd\nchrI\t1\n")


def test_float_mutation_summary_does_not_crash(tmp_path):
    path = tmp_path / "mut.tsv"
    _mutations(path)
    metrics = load_metric_table(str(path))
    assert metrics["forward_fraction"] == 0.625


def test_composition_depends_on_the_counts(tmp_path):
    fasta = tmp_path / "genome.fa"
    pathways = tmp_path / "pathways.yaml"
    mutations = tmp_path / "mut.tsv"
    _fasta(fasta)
    _pathways(pathways)
    _mutations(mutations)
    balanced = tmp_path / "balanced.tsv"
    skewed = tmp_path / "skewed.tsv"
    _signal(
        balanced,
        [("chrI", 0, "+", 10, "A"), ("chrI", 1, "+", 10, "C"), ("chrI", 2, "+", 10, "G"), ("chrI", 3, "+", 10, "T")],
    )
    _signal(skewed, [("chrI", 0, "+", 40, "A")])
    even = run(str(balanced), str(mutations), str(pathways), str(fasta), None, False)
    odd = run(str(skewed), str(mutations), str(pathways), str(fasta), None, False)
    even_stat = even["signal_tests"][0]["statistic"]
    odd_stat = odd["signal_tests"][0]["statistic"]
    assert even["signal_tests"][0]["id"] == "rNMP_base_composition"
    assert even_stat != odd_stat
    assert odd_stat > even_stat
    assert all(row["enrichment_score"] is None for row in even["pathways"])
    assert even["pathway_ora"] == "disabled"


def test_hypergeometric_increases_when_overlap_increases():
    low = hypergeometric_pvalue(1, 10, 10, 100)
    high = hypergeometric_pvalue(8, 10, 10, 100)
    assert high < low
    assert 0.0 <= high <= 1.0


def test_ora_uses_gene_overlap_not_panel_size(tmp_path):
    fasta = tmp_path / "genome.fa"
    pathways = tmp_path / "pathways.yaml"
    mutations = tmp_path / "mut.tsv"
    signal = tmp_path / "signal.tsv"
    bed = tmp_path / "genes.bed"
    _fasta(fasta)
    _pathways(pathways)
    _mutations(mutations)
    _signal(signal, [("chrI", 10, "+", 3, "A"), ("chrI", 500, "+", 1, "C")])
    bed.write_text(
        "chrI\t0\t20\tGENE1\tpathway-a\n"
        "chrI\t100\t120\tGENE2\tpathway-a\n"
        "chrI\t200\t220\tGENE3\tpathway-b\n"
        "chrI\t480\t520\tGENE4\tpathway-b\n"
    )
    payload = run(str(signal), str(mutations), str(pathways), str(fasta), str(bed), True)
    by_id = {row["pathway_id"]: row for row in payload["pathways"]}
    assert payload["pathway_ora"] == "computed"
    assert by_id["pathway-a"]["overlap"] == 1
    assert by_id["pathway-b"]["overlap"] == 1
    assert by_id["pathway-a"]["enrichment_score"] != 0.1
