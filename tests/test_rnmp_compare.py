"""Wild-type versus rnh201 comparison of direct rNMP maps."""

import csv
import json
from pathlib import Path

from rnmp_compare import compare_samples, write_html

ROOT = Path(__file__).resolve().parents[1]


def _signature(accession, genotype, ends, positions, per_mb, forward, reverse, bases):
    total = forward + reverse
    return {
        "accession": accession,
        "genotype": genotype,
        "technique": "ribose-seq",
        "assay_type": "ribose_seq",
        "total_rnmp_5prime_ends": ends,
        "unique_rnmp_positions": positions,
        "rnmp_calls_per_mb": per_mb,
        "forward_ends": forward,
        "reverse_ends": reverse,
        "forward_fraction": round(forward / total, 4),
        "avg_strand_bias": 0.5,
        "base_counts": bases,
        "polymerase_assignment": "experimental_not_assigned",
    }


def test_compare_reports_counts_strand_composition_and_enrichment():
    wild = _signature("SRR11364941", "wild_type", 10, 4, 1.0, 7, 3, {"A": 6, "C": 4, "G": 0, "T": 0})
    mutant = _signature("SRR11364929", "rnh201", 40, 16, 4.0, 10, 30, {"A": 10, "C": 30, "G": 0, "T": 0})
    payload = compare_samples(
        [wild, mutant],
        {
            "SRR11364941": {
                "accession": "SRR11364941",
                "signal_tests": [{"id": "strand_balance", "status": "computed", "statistic": 0.7, "pvalue": 0.3}],
            }
        },
    )
    assert payload["status"] == "compared"
    assert payload["position_counts"]["wild_type_total_ends"] == 10
    assert payload["position_counts"]["rnh201_unique_positions"] == 16
    assert payload["enrichment"]["total_ends_rnh201_over_wt"] == 4.0
    assert payload["enrichment"]["unique_positions_rnh201_over_wt"] == 4.0
    assert payload["enrichment"]["calls_per_mb_rnh201_over_wt"] == 4.0
    assert payload["nucleotide_composition"]["status"] == "computed"
    assert payload["strand_bias"]["status"] == "computed"
    assert payload["strand_bias"]["forward_fraction_difference_rnh201_minus_wt"] == round(10 / 40 - 7 / 10, 4)
    assert payload["wild_type"]["signal_tests"]["strand_balance"]["pvalue"] == 0.3
    assert payload["polymerase_assignment"] == "experimental_not_assigned"
    assert "qc" not in payload


def test_compare_skips_without_both_genotypes():
    only = _signature("SRR1575904", "", 5, 5, 1, 1, 4, {"A": 5, "C": 0, "G": 0, "T": 0})
    payload = compare_samples([only])
    assert payload["status"] == "skipped"


def test_comparison_page_keeps_v1_scope(tmp_path):
    wild = _signature("SRR11364941", "wild_type", 10, 4, 1.0, 7, 3, {"A": 6, "C": 4, "G": 0, "T": 0})
    mutant = _signature("SRR11364929", "rnh201", 40, 16, 4.0, 10, 30, {"A": 10, "C": 30, "G": 0, "T": 0})
    payload = compare_samples([wild, mutant])
    page = tmp_path / "comparison_report.html"
    write_html(payload, page)
    text = page.read_text()
    before, future = text.split('id="future-work"', 1)
    assert "SRR11364941" in before
    assert "SRR11364929" in before
    assert "hydroxyurea" not in before.lower()
    assert "T>C" not in text
    assert "experimental_not_assigned" in before
    assert "30847483" in future
    assert "25664722" in future
    json.dumps(payload)


def test_v1_sample_sheet_is_the_wt_rnh201_pair():
    with (ROOT / "samples_test.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert [row["accession"] for row in rows] == ["SRR11364941", "SRR11364929"]
    assert [row["genotype"] for row in rows] == ["wild_type", "rnh201"]
    assert {row["technique"] for row in rows} == {"ribose-seq"}
    assert {row["trim_5prime"] for row in rows} == {"12"}
    assert all("32415081" in row["notes"] for row in rows)
    assert all(row["adapter_3prime"] == "AGTTGCGACACGGATCTATCA" for row in rows)


def test_browser_presents_v1_not_an_hu_result():
    text = (ROOT / "docker/browser/public/index.html").read_text()
    before, future = text.split('id="future-work"', 1)
    assert "HU-rNMP" not in text
    assert "hydroxyurea" not in before.lower()
    assert "wild-type" in before.lower() or "wild type" in before.lower()
    assert "RNase H2" in before
    assert "30847483" in future
    assert "25664722" in future
    assert "Pol ε" not in text
    assert "Pol alpha" not in text
