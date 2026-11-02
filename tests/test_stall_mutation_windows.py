"""Window intersection: VCF 1-based positions vs 0-based HydEn bins, plus promoters."""

from stall_mutation_windows import (
    C0,
    HU,
    annotate_sites,
    bin_start,
    build_windows,
    end_proximal,
    index_genes,
    promoter_span,
    site_transcription,
    wgs_pos_to_zero,
)


def _gene(**kwargs):
    base = {
        "chrom": "chrI",
        "start": 200,
        "end": 400,
        "strand": "+",
        "systematic": "YAL000W",
        "symbol": "GENE",
        "biotype": "protein_coding",
    }
    base.update(kwargs)
    return base


def test_end_proximal_flags_chromosome_ends_and_mt():
    lengths = {"chrI": 230218}
    assert end_proximal("chrI", 0, 500, lengths) == 1
    assert end_proximal("chrI", 100000, 100500, lengths) == 0
    assert end_proximal("chrI", 220000, 220500, lengths) == 1
    assert end_proximal("chrmt", 29000, 29500, lengths) == 1


def test_vcf_boundary_stays_in_lower_bin():
    assert wgs_pos_to_zero(500) == 499
    assert bin_start(499) == 0
    assert bin_start(500) == 500
    assert bin_start(wgs_pos_to_zero(501)) == 500


def test_promoter_is_strand_aware():
    plus = _gene(strand="+", start=200, end=400)
    minus = _gene(strand="-", start=1000, end=1500, systematic="YAL000C", symbol="REV")
    assert promoter_span(plus) == (-300, 200)
    assert promoter_span(minus) == (1500, 2000)
    klass, bodies, promoters = site_transcription(99, [plus])
    assert klass == "promoter"
    assert bodies == []
    assert promoters == [plus]
    klass, bodies, promoters = site_transcription(249, [plus])
    assert klass == "gene_body"
    assert bodies == [plus]
    klass, _, promoters = site_transcription(1600, [minus])
    assert klass == "promoter"
    assert promoters == [minus]


def test_continued_stall_requires_both_libraries_and_transcription():
    genes = index_genes([_gene()])
    nt = {("chrI", 0): {"chrom": "chrI", "start": 0, "end": 500, "total": 100, "fold": 4.0}}
    hii = {
        ("chrI", 0): {"chrom": "chrI", "start": 0, "end": 500, "total": 80, "fold": 3.0},
        ("chrI", 500): {"chrom": "chrI", "start": 500, "end": 1000, "total": 50, "fold": 2.5},
    }
    sites = [
        {"run": HU, "chrom": "chrI", "pos": 100, "pos0": 99, "ref": "T", "alt": "C"},
        {"run": C0, "chrom": "chrI", "pos": 100, "pos0": 99, "ref": "T", "alt": "C"},
        {"run": HU, "chrom": "chrI", "pos": 250, "pos0": 249, "ref": "A", "alt": "G"},
        {"run": HU, "chrom": "chrI", "pos": 600, "pos0": 599, "ref": "T", "alt": "C"},
        {"run": HU, "chrom": "chrI", "pos": 2000, "pos0": 1999, "ref": "T", "alt": "C"},
        {"run": C0, "chrom": "chrI", "pos": 500, "pos0": 499, "ref": "A", "alt": "G"},
    ]
    shared = {("chrI", 100, "T", "C")}
    annotated = annotate_sites(sites, nt, hii, genes, shared)
    by_pos = {(row["run"], row["pos"]): row for row in annotated}

    boundary = by_pos[(C0, 500)]
    assert boundary["bin_start"] == 0
    assert boundary["hotspot_class"] == "continued"
    assert boundary["tx_class"] == "intergenic"
    assert boundary["continued_stall"] == 0

    promoter = by_pos[(HU, 100)]
    assert promoter["tx_class"] == "promoter"
    assert promoter["continued_stall"] == 1
    assert promoter["shared_site"] == 1

    hii_only = by_pos[(HU, 600)]
    assert hii_only["hotspot_class"] == "hii_only"
    assert hii_only["continued_stall"] == 0

    outside = by_pos[(HU, 2000)]
    assert outside["hotspot_class"] == "none"
    assert outside["tx_class"] == "intergenic"

    windows = build_windows(annotated)
    continued = [row for row in windows if row["start"] == 0]
    assert len(continued) == 1
    assert continued[0]["continued_stall"] == 1
    assert continued[0]["hotspot_class"] == "continued"
    assert continued[0]["n_hu"] == 2
    assert continued[0]["n_c0"] == 2
    assert continued[0]["n_shared"] == 1
    hii_windows = [row for row in windows if row["start"] == 500]
    assert hii_windows[0]["continued_stall"] == 0
    assert hii_windows[0]["hotspot_class"] == "hii_only"
    assert all(row["start"] != 1500 for row in windows)
