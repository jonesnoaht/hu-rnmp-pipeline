"""SGD intersection: overlap shares, chrM alias, and the HII/NT ratio rule."""

from gene_excision_rank import (
    _by_chrom,
    classify_excision,
    genes_overlapping,
    hotspot_shares,
    norm_chrom,
    parse_gff_genes,
)


def _gene(**kwargs):
    base = {
        "chrom": "chrI",
        "start": 400,
        "end": 600,
        "strand": "+",
        "systematic": "YAL000W",
        "symbol": "GENE",
        "biotype": "gene",
        "dubious": False,
        "te": False,
        "rrna": False,
        "mitochondrial": False,
    }
    base.update(kwargs)
    return base


def test_chrmt_aliases_to_pipeline_chrM():
    assert norm_chrom("chrmt") == "chrM"
    assert norm_chrom("chrXIV") == "chrXIV"


def test_nt_zero_is_not_a_ratio():
    klass, ratio = classify_excision(0, 2612, 0, 40)
    assert klass == "HII_only"
    assert ratio is None


def test_supported_ratio_is_unique_positions_not_depth():
    klass, ratio = classify_excision(1000, 800, 40, 100)
    assert klass == "HII_enriched"
    assert ratio == 2.5


def test_end_pileup_does_not_set_the_class():
    klass, ratio = classify_excision(1000, 10000, 50, 55)
    assert klass == "shared"
    assert ratio == 1.1


def test_below_floor_is_sparse_even_if_the_site_ratio_is_large():
    klass, ratio = classify_excision(10, 80, 4, 12)
    assert klass == "sparse"
    assert ratio == 3.0


def test_hotspot_mass_follows_overlap_not_an_equal_split():
    gene_a = _gene(start=0, end=200, systematic="A")
    gene_b = _gene(start=200, end=250, systematic="B")
    shares = {
        gene["systematic"]: share
        for gene, share in hotspot_shares(0, 500, 100, [gene_a, gene_b])
    }
    assert shares["A"] == 40
    assert shares["B"] == 10
    assert shares["A"] != 50


def test_point_hits_transcript_extended_span_and_chrM():
    nuclear = _gene(start=100, end=120, systematic="YAL001W")
    mito = _gene(chrom="chrM", start=10, end=30, systematic="Q0045", symbol="COX1", mitochondrial=True)
    indexed = _by_chrom([nuclear, mito])
    assert [g["systematic"] for g in genes_overlapping(indexed["chrI"], 110, 111)] == ["YAL001W"]
    assert genes_overlapping(indexed["chrI"], 120, 121) == []
    assert [g["systematic"] for g in genes_overlapping(indexed["chrM"], 10, 11)] == ["Q0045"]


def test_gff_parser_uses_gene_symbol_and_zero_based_span(tmp_path):
    gff = tmp_path / "mini.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chrVI\tSGD\tchromosome\t1\t270161\t.\t.\t.\tID=chrVI;Name=chrVI\n"
        "chrmt\tSGD\tchromosome\t1\t85779\t.\t.\t.\tID=chrmt;Name=chrmt\n"
        "chrVI\tSGD\tgene\t52971\t54910\t.\t-\t.\tID=YFL039C;Name=YFL039C;gene=ACT1;Note=Actin\n"
        "chrVI\tSGD\tCDS\t52971\t54465\t.\t-\t.\tID=CDS:YFL039C;Parent=YFL039C\n"
        "chrmt\tSGD\tgene\t100\t200\t.\t+\t.\tID=Q0045;Name=Q0045;gene=COX1;Note=cytochrome\n"
        "chrI\tSGD\tgene\t335\t649\t.\t-\t.\tID=YAL069W;Name=YAL069W;Note=Dubious%20open%20reading%20frame\n"
    )
    genes, lengths = parse_gff_genes(gff)
    by_id = {g["systematic"]: g for g in genes}
    assert set(by_id) == {"YFL039C", "Q0045", "YAL069W"}
    assert by_id["YFL039C"]["symbol"] == "ACT1"
    assert by_id["YFL039C"]["start"] == 52970
    assert by_id["YFL039C"]["end"] == 54910
    assert by_id["Q0045"]["chrom"] == "chrM"
    assert by_id["YAL069W"]["dubious"] is True
    assert lengths["chrM"] == 85779
    assert lengths["chrVI"] == 270161
