#!/usr/bin/env python3
"""Signal statistics, plus optional gene-set over-representation.

The previous pathway score was ``len(gene_panel) * 0.1``. It did not read the
rNMP table. That score is gone.

What runs by default is a data-dependent goodness-of-fit:

- chi-square of strand-aware rNMP base counts against the reference
  mononucleotide composition (A paired with T, C paired with G, because a
  minus-strand call is the complement of the reference base);
- a two-sided binomial test of plus versus minus call counts.

Neither test is a pathway enrichment. Ranking the human gene panels against
yeast rNMP positions is not run unless ``--enable-pathway-ora`` and a
``--genes-bed`` file are both set. The BED columns are
``chrom start end gene pathway_id`` (0-based, half-open). The test is a
one-sided hypergeometric over-representation of genes overlapped by at least
one rNMP call.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter

import yaml
from scipy import stats

BASES = ("A", "C", "G", "T")


def load_pathways(pathways_yaml: str) -> dict:
    with open(pathways_yaml) as handle:
        data = yaml.safe_load(handle)
    return {row["id"]: row for row in data.get("pathways", [])}


def load_direct_signal(signal_tsv: str) -> list[dict]:
    """Load HydEn/ribose signal rows. Returns [] for a non-matching header."""
    rows: list[dict] = []
    with open(signal_tsv) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: pos for pos, name in enumerate(header)}
        if "strand" not in index or "count" not in index:
            return []
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= index["count"]:
                continue
            try:
                count = int(parts[index["count"]])
                pos = int(parts[index["pos"]]) if "pos" in index else None
            except ValueError:
                continue
            base = parts[index["base"]] if "base" in index and index["base"] < len(parts) else "NA"
            rows.append(
                {
                    "chrom": parts[index["chrom"]] if "chrom" in index else "",
                    "pos": pos,
                    "strand": parts[index["strand"]],
                    "count": count,
                    "base": base.upper(),
                }
            )
    return rows


def load_metric_table(mutations_tsv: str) -> dict:
    """Read a two-column metric table. Floats are kept; a second table is ignored.

    The old HydEn summary wrote ``forward_fraction\\t0.6250`` and then a second
    header. ``int()`` on that value crashed the enrichment stage.
    """
    metrics: dict = {}
    with open(mutations_tsv) as handle:
        handle.readline()
        for line in handle:
            if not line.strip():
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "chrom":
                break
            key, raw = parts[0], parts[1]
            try:
                metrics[key] = int(raw)
            except ValueError:
                try:
                    metrics[key] = float(raw)
                except ValueError:
                    metrics[key] = raw
    return metrics


def reference_base_frequencies(fasta_path: str) -> dict[str, float]:
    counts: Counter[str] = Counter()
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            counts.update(base for base in line.strip().upper() if base in BASES)
    total = sum(counts[base] for base in BASES)
    if total == 0:
        raise ValueError(f"no A/C/G/T bases in {fasta_path}")
    paired = {
        "A": counts["A"] + counts["T"],
        "T": counts["A"] + counts["T"],
        "C": counts["C"] + counts["G"],
        "G": counts["C"] + counts["G"],
    }
    denom = float(sum(paired.values()))
    return {base: paired[base] / denom for base in BASES}


def base_composition_test(rows: list[dict], frequencies: dict[str, float]) -> dict:
    observed = {base: 0 for base in BASES}
    for row in rows:
        if row["base"] in observed:
            observed[row["base"]] += row["count"]
    total = sum(observed.values())
    if total == 0:
        return {
            "id": "rNMP_base_composition",
            "status": "not_computed",
            "statistic": None,
            "pvalue": None,
            "observed": observed,
            "expected_frequency": frequencies,
            "note": "No A/C/G/T calls were available to test.",
        }
    obs = [observed[base] for base in BASES]
    exp = [frequencies[base] * total for base in BASES]
    statistic, pvalue = stats.chisquare(obs, exp)
    return {
        "id": "rNMP_base_composition",
        "status": "computed",
        "test": "chisquare",
        "statistic": round(float(statistic), 6),
        "pvalue": float(pvalue),
        "observed": observed,
        "expected_frequency": {base: round(frequencies[base], 6) for base in BASES},
        "note": (
            "Goodness-of-fit of strand-aware rNMP base counts versus the "
            "reference mononucleotide composition. This is not a pathway enrichment "
            "and is not corrected for sequence context or mappability."
        ),
    }


def strand_balance_test(rows: list[dict]) -> dict:
    plus = sum(row["count"] for row in rows if row["strand"] == "+")
    minus = sum(row["count"] for row in rows if row["strand"] == "-")
    total = plus + minus
    if total == 0:
        return {
            "id": "strand_balance",
            "status": "not_computed",
            "statistic": None,
            "pvalue": None,
            "observed": {"+": plus, "-": minus},
            "note": "No stranded calls were available to test.",
        }
    pvalue = float(stats.binomtest(plus, total, 0.5, alternative="two-sided").pvalue)
    return {
        "id": "strand_balance",
        "status": "computed",
        "test": "binomial",
        "statistic": round(plus / total, 6),
        "pvalue": pvalue,
        "observed": {"+": plus, "-": minus},
        "expected": 0.5,
        "note": (
            "Two-sided binomial test of plus versus minus rNMP calls. "
            "A low p-value is strand imbalance, not a polymerase assignment."
        ),
    }


def hypergeometric_pvalue(overlap: int, pathway_size: int, n_hits: int, universe: int) -> float:
    """One-sided P(X >= overlap) for a Hypergeometric draw."""
    if min(overlap, pathway_size, n_hits, universe) < 0:
        raise ValueError("hypergeometric counts must be non-negative")
    if universe == 0 or pathway_size == 0 or n_hits == 0:
        return 1.0
    if overlap > min(pathway_size, n_hits):
        return 0.0
    return float(stats.hypergeom.sf(overlap - 1, universe, pathway_size, n_hits))


def overlap_genes(rows: list[dict], genes_bed: str) -> tuple[set[str], dict[str, set[str]]]:
    """Return hit genes and pathway -> genes from a 5-column BED."""
    genes: list[tuple[str, int, int, str, str]] = []
    with open(genes_bed) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "chrom")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            genes.append((parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]))
    by_chrom: dict[str, list[tuple[int, int, str, str]]] = {}
    pathway_genes: dict[str, set[str]] = {}
    for chrom, start, end, gene, pathway in genes:
        by_chrom.setdefault(chrom, []).append((start, end, gene, pathway))
        pathway_genes.setdefault(pathway, set()).add(gene)
    hits: set[str] = set()
    for row in rows:
        if row["pos"] is None or row["count"] <= 0:
            continue
        for start, end, gene, _pathway in by_chrom.get(row["chrom"], []):
            if start <= row["pos"] < end:
                hits.add(gene)
    return hits, pathway_genes


def pathway_rows_disabled(pathways: dict) -> list[dict]:
    note = (
        "Over-representation was not run. Supply --genes-bed and "
        "--enable-pathway-ora to test gene overlap with a hypergeometric test. "
        "Human gene-panel size is not a score for these positions."
    )
    rows = []
    for pid, meta in pathways.items():
        rows.append(
            {
                "pathway_id": pid,
                "pathway_name": meta.get("name", pid),
                "status": "not_computed",
                "enrichment_score": None,
                "pvalue": None,
                "overlap": None,
                "pathway_size": None,
                "note": note,
            }
        )
    return rows


def pathway_ora(rows: list[dict], genes_bed: str) -> list[dict]:
    hits, pathway_genes = overlap_genes(rows, genes_bed)
    universe = {gene for genes in pathway_genes.values() for gene in genes}
    n_universe = len(universe)
    n_hits = len(hits & universe)
    results = []
    for pathway, genes in sorted(pathway_genes.items()):
        overlap = len(hits & genes)
        size = len(genes)
        expected = n_hits * (size / n_universe) if n_universe else 0.0
        pvalue = hypergeometric_pvalue(overlap, size, n_hits, n_universe)
        fold = (overlap / expected) if expected else None
        results.append(
            {
                "pathway_id": pathway,
                "pathway_name": pathway,
                "status": "computed",
                "enrichment_score": None if fold is None else round(fold, 6),
                "pvalue": pvalue,
                "overlap": overlap,
                "pathway_size": size,
                "n_hits": n_hits,
                "universe": n_universe,
                "note": (
                    "Hypergeometric over-representation of genes overlapped by an rNMP call. "
                    "enrichment_score is the fold over the expected overlap, not a panel-size weight."
                ),
            }
        )
    results.sort(key=lambda row: (row["pvalue"], -(row["overlap"] or 0)))
    return results


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, dict):
        return ";".join(f"{key}:{val}" for key, val in value.items())
    return str(value).replace("\t", " ").replace("\n", " ")


def write_outputs(payload: dict, output_tsv: str, output_json: str) -> None:
    with open(output_tsv, "w") as handle:
        handle.write(
            "record_type\tid\tname\tstatus\tstatistic\tpvalue\tobserved\texpected\tnote\n"
        )
        lines = []
        for test in payload["signal_tests"]:
            lines.append(
                "\t".join(
                    [
                        "signal_test",
                        test["id"],
                        test.get("test", ""),
                        test["status"],
                        _fmt(test.get("statistic")),
                        _fmt(test.get("pvalue")),
                        _fmt(test.get("observed")),
                        _fmt(test.get("expected_frequency", test.get("expected"))),
                        _fmt(test.get("note")),
                    ]
                )
                + "\n"
            )
        for row in payload["pathways"]:
            lines.append(
                "\t".join(
                    [
                        "pathway",
                        row["pathway_id"],
                        _fmt(row.get("pathway_name")),
                        row["status"],
                        _fmt(row.get("enrichment_score")),
                        _fmt(row.get("pvalue")),
                        _fmt(row.get("overlap")),
                        _fmt(row.get("pathway_size")),
                        _fmt(row.get("note")),
                    ]
                )
                + "\n"
            )
        handle.writelines(lines)
    with open(output_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def run(
    signal_tsv: str,
    mutations_tsv: str,
    pathways_yaml: str,
    fasta_path: str | None,
    genes_bed: str | None,
    enable_pathway_ora: bool,
) -> dict:
    pathways = load_pathways(pathways_yaml)
    signal_rows = load_direct_signal(signal_tsv)
    metrics = load_metric_table(mutations_tsv) if mutations_tsv else {}
    tests = [strand_balance_test(signal_rows)]
    if fasta_path and signal_rows:
        tests.insert(0, base_composition_test(signal_rows, reference_base_frequencies(fasta_path)))
    elif fasta_path:
        tests.insert(
            0,
            {
                "id": "rNMP_base_composition",
                "status": "not_computed",
                "statistic": None,
                "pvalue": None,
                "note": "Signal table had no direct-assay rows.",
            },
        )
    if enable_pathway_ora and genes_bed:
        pathway_rows = pathway_ora(signal_rows, genes_bed)
        ora_status = "computed"
    else:
        pathway_rows = pathway_rows_disabled(pathways)
        ora_status = "disabled"
    return {
        "pathway_ora": ora_status,
        "signal_tests": tests,
        "pathways": pathway_rows,
        "total_signal_calls": sum(row["count"] for row in signal_rows),
        "input_metrics": metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="rNMP signal tests and optional pathway ORA")
    parser.add_argument("--signal", required=True)
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--pathways", required=True)
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--genes-bed", default=None)
    parser.add_argument("--enable-pathway-ora", action="store_true")
    parser.add_argument("--model-packs", default=None, help="Accepted and unused")
    parser.add_argument("--output", required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    payload = run(
        args.signal,
        args.mutations,
        args.pathways,
        args.fasta,
        args.genes_bed,
        args.enable_pathway_ora,
    )
    write_outputs(payload, args.output, args.json)
    computed = [row["id"] for row in payload["signal_tests"] if row["status"] == "computed"]
    print(
        f"Signal tests computed: {', '.join(computed) or 'none'}; "
        f"pathway ORA {payload['pathway_ora']} → {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
