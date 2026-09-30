"""ENA filereport parsing and reproducible subsampling. No network."""

import gzip

import pytest
from fetch_ena import (
    fastq_entries,
    ftp_to_https,
    parse_filereport,
    select_run,
    subsample_fastq,
)

REPORT = (
    "run_accession\tsample_title\tlibrary_layout\tfastq_ftp\tfastq_md5\tfastq_bytes\tread_count\n"
    "SRR1575904\tpol2-4 rnh201 (KK-107)\tSINGLE\t"
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1575904/SRR1575904.fastq.gz\t"
    "0d3bdfeda49f790f49edc83bb1a47050\t39337368\t1068754\n"
)


def test_parse_single_end_report():
    rows = parse_filereport(REPORT)
    assert rows[0]["run_accession"] == "SRR1575904"
    assert rows[0]["library_layout"] == "SINGLE"
    url, digest = fastq_entries(rows[0])[0]
    assert url == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1575904/SRR1575904.fastq.gz"
    assert digest == "0d3bdfeda49f790f49edc83bb1a47050"


def test_ftp_url_schemes():
    assert ftp_to_https("ftp://ftp.sra.ebi.ac.uk/vol1/a.fastq.gz") == "https://ftp.sra.ebi.ac.uk/vol1/a.fastq.gz"
    assert ftp_to_https("http://ftp.sra.ebi.ac.uk/a") == "https://ftp.sra.ebi.ac.uk/a"


def test_study_is_not_expanded_unless_asked():
    rows = [
        {"run_accession": "SRR1", "fastq_bytes": "50"},
        {"run_accession": "SRR2", "fastq_bytes": "10"},
    ]
    with pytest.raises(ValueError, match="Refusing"):
        select_run(rows, "SRP182841", max_runs=0)
    assert select_run(rows, "SRP182841", max_runs=1)["run_accession"] == "SRR2"


def test_gse_is_rejected():
    with pytest.raises(ValueError, match="GEO series"):
        select_run([], "GSE61464", max_runs=0)


def test_subsample_is_deterministic(tmp_path):
    src = tmp_path / "in.fastq.gz"
    with gzip.open(src, "wt") as handle:
        for index in range(20):
            handle.write(f"@r{index}\n{'A' * 8}{index:02d}\n+\n{'I' * 10}\n")
    first = tmp_path / "a.fastq.gz"
    second = tmp_path / "b.fastq.gz"
    assert subsample_fastq(src, first, 5, 1) == 5
    subsample_fastq(src, second, 5, 1)

    def _text(path):
        with gzip.open(path, "rt") as handle:
            return handle.read()

    assert _text(first) == _text(second)
    other = tmp_path / "c.fastq.gz"
    subsample_fastq(src, other, 5, 2)
    assert _text(first) != _text(other)
