#!/usr/bin/env python3
"""Compare wild-type and rnh201 direct rNMP maps.

The contrast is position counts, strand bias, nucleotide composition, and
enrichment of the RNase H2 mutant over wild type. It is a description of the
libraries that were passed in, not a multi-replicate significance claim.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path

from build_report import V1_SCOPE, future_work_section
from scipy import stats

WILD_TYPE = "wild_type"
MUTANT = "rnh201"
BASES = ("A", "C", "G", "T")


def _load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise TypeError(f"{path} is not a JSON object")
    return payload


def _json_files(path: Path) -> list[Path]:
    if path.is_dir():
        return sorted(item for item in path.glob("*.json") if item.is_file())
    if path.is_file():
        return [path]
    return []


def _esc(value) -> str:
    return html.escape("" if value is None else str(value))


def _fraction(counts: dict, key: str, total: int) -> float | None:
    if total <= 0:
        return None
    return counts.get(key, 0) / total


def _ratio(numerator, denominator):
    if numerator is None or denominator in (None, 0):
        return None
    return round(numerator / denominator, 4)


def _sample_view(signature: dict, hits: dict | None) -> dict:
    bases = {base: int(signature.get("base_counts", {}).get(base, 0) or 0) for base in BASES}
    total_bases = sum(bases.values())
    fractions = {}
    for base in BASES:
        frac = _fraction(bases, base, total_bases)
        fractions[base] = None if frac is None else round(frac, 4)
    forward = int(signature.get("forward_ends") or 0)
    reverse = int(signature.get("reverse_ends") or 0)
    total_ends = int(signature.get("total_rnmp_5prime_ends") or (forward + reverse))
    tests = {}
    for test in (hits or {}).get("signal_tests") or []:
        if test.get("id"):
            tests[test["id"]] = {
                "status": test.get("status"),
                "statistic": test.get("statistic"),
                "pvalue": test.get("pvalue"),
            }
    return {
        "accession": signature.get("accession") or "",
        "genotype": signature.get("genotype") or "",
        "technique": signature.get("technique") or "",
        "assay_type": signature.get("assay_type") or "",
        "total_rnmp_5prime_ends": total_ends,
        "unique_rnmp_positions": signature.get("unique_rnmp_positions"),
        "rnmp_calls_per_mb": signature.get("rnmp_calls_per_mb"),
        "forward_ends": forward,
        "reverse_ends": reverse,
        "forward_fraction": signature.get("forward_fraction"),
        "avg_strand_bias": signature.get("avg_strand_bias"),
        "base_counts": bases,
        "base_fractions": fractions,
        "signal_tests": tests,
    }


def _composition_contrast(wild: dict, mutant: dict) -> dict:
    bases = [base for base in BASES if wild["base_counts"][base] + mutant["base_counts"][base] > 0]
    wild_total = sum(wild["base_counts"].values())
    mutant_total = sum(mutant["base_counts"].values())
    if len(bases) < 2 or wild_total == 0 or mutant_total == 0:
        return {
            "status": "not_computed",
            "note": "Need counts on at least two bases in both genotypes.",
        }
    table = [
        [wild["base_counts"][base] for base in bases],
        [mutant["base_counts"][base] for base in bases],
    ]
    statistic, pvalue, dof, _expected = stats.chi2_contingency(table)
    return {
        "status": "computed",
        "test": "chi2_contingency",
        "bases": bases,
        "statistic": round(float(statistic), 6),
        "dof": int(dof),
        "pvalue": float(pvalue),
        "note": (
            "Contingency test of incorporated-base counts between genotypes. "
            "Not corrected for sequence context, mappability, or library depth beyond the raw counts."
        ),
    }


def _strand_contrast(wild: dict, mutant: dict) -> dict:
    table = [
        [wild["forward_ends"], wild["reverse_ends"]],
        [mutant["forward_ends"], mutant["reverse_ends"]],
    ]
    if any(sum(row) == 0 for row in table) or any(table[0][i] + table[1][i] == 0 for i in (0, 1)):
        return {"status": "not_computed", "note": "Both genotypes need plus and minus counts."}
    statistic, pvalue = stats.fisher_exact(table)
    wild_fraction = wild["forward_fraction"]
    mutant_fraction = mutant["forward_fraction"]
    delta = None
    if wild_fraction is not None and mutant_fraction is not None:
        delta = round(mutant_fraction - wild_fraction, 4)
    return {
        "status": "computed",
        "test": "fisher_exact",
        "statistic": round(float(statistic), 6),
        "pvalue": float(pvalue),
        "forward_fraction_difference_rnh201_minus_wt": delta,
        "note": (
            "Fisher exact test of strand by genotype. A low p-value is a difference in "
            "strand balance. v1 does not assign replicative polymerases."
        ),
    }


def compare_samples(signatures: list[dict], hits_by_accession: dict[str, dict] | None = None) -> dict:
    """Build the wild-type versus rnh201 comparison, or a skipped payload."""
    hits_by_accession = hits_by_accession or {}
    views = []
    for signature in signatures:
        accession = signature.get("accession") or ""
        views.append(_sample_view(signature, hits_by_accession.get(accession)))
    wild = [row for row in views if row["genotype"] == WILD_TYPE]
    mutant = [row for row in views if row["genotype"] == MUTANT]
    if len(wild) != 1 or len(mutant) != 1:
        return {
            "status": "skipped",
            "reason": (
                f"Need exactly one {WILD_TYPE} and one {MUTANT} signature; "
                f"found {len(wild)} and {len(mutant)}."
            ),
            "samples": [
                {"accession": row["accession"], "genotype": row["genotype"]} for row in views
            ],
        }
    left, right = wild[0], mutant[0]
    enrichment = {
        "total_ends_rnh201_over_wt": _ratio(right["total_rnmp_5prime_ends"], left["total_rnmp_5prime_ends"]),
        "unique_positions_rnh201_over_wt": _ratio(right["unique_rnmp_positions"], left["unique_rnmp_positions"]),
        "calls_per_mb_rnh201_over_wt": _ratio(right["rnmp_calls_per_mb"], left["rnmp_calls_per_mb"]),
        "note": (
            "Ratio of rnh201 over wild type. One library per genotype, using the counts "
            "this run produced (including a subsample when the test profile set one). "
            "Not a replicate-aware enrichment model."
        ),
    }
    return {
        "status": "compared",
        "comparison": "wild_type_vs_rnh201",
        "scope": V1_SCOPE,
        "metrics": [
            "rnmp_position_counts",
            "strand_bias",
            "nucleotide_composition",
            "enrichment",
        ],
        "wild_type": left,
        "rnh201": right,
        "position_counts": {
            "wild_type_total_ends": left["total_rnmp_5prime_ends"],
            "rnh201_total_ends": right["total_rnmp_5prime_ends"],
            "wild_type_unique_positions": left["unique_rnmp_positions"],
            "rnh201_unique_positions": right["unique_rnmp_positions"],
            "wild_type_calls_per_mb": left["rnmp_calls_per_mb"],
            "rnh201_calls_per_mb": right["rnmp_calls_per_mb"],
        },
        "strand_bias": _strand_contrast(left, right),
        "nucleotide_composition": _composition_contrast(left, right),
        "enrichment": enrichment,
        "mutation_signature_note": (
            "The 2–5 bp repeat-deletion counter is the headline mutation signature for "
            "variant calls. This comparison is a direct rNMP map, not a deletion callset."
        ),
        "polymerase_assignment": "experimental_not_assigned",
    }


def _index_hits(payloads: list[dict]) -> dict[str, dict]:
    indexed = {}
    for payload in payloads:
        accession = payload.get("accession") or ""
        if accession:
            indexed[accession] = payload
    return indexed


def write_table(payload: dict, path: Path) -> None:
    lines = ["section\tmetric\twild_type\trnh201\tvalue\n"]
    if payload.get("status") != "compared":
        lines.append(f"status\tskipped\t\t\t{payload.get('reason', '')}\n")
        path.write_text("".join(lines))
        return
    counts = payload["position_counts"]
    lines.extend(
        [
            f"position_counts\ttotal_ends\t{counts['wild_type_total_ends']}\t{counts['rnh201_total_ends']}\t\n",
            (
                "position_counts\tunique_positions\t"
                f"{counts['wild_type_unique_positions']}\t{counts['rnh201_unique_positions']}\t\n"
            ),
            (
                "position_counts\tcalls_per_mb\t"
                f"{counts['wild_type_calls_per_mb']}\t{counts['rnh201_calls_per_mb']}\t\n"
            ),
            (
                "enrichment\ttotal_ends_rnh201_over_wt\t\t\t"
                f"{payload['enrichment']['total_ends_rnh201_over_wt']}\n"
            ),
            (
                "enrichment\tunique_positions_rnh201_over_wt\t\t\t"
                f"{payload['enrichment']['unique_positions_rnh201_over_wt']}\n"
            ),
            (
                "enrichment\tcalls_per_mb_rnh201_over_wt\t\t\t"
                f"{payload['enrichment']['calls_per_mb_rnh201_over_wt']}\n"
            ),
        ]
    )
    strand = payload["strand_bias"]
    lines.append(
        "strand_bias\tforward_fraction_difference_rnh201_minus_wt\t"
        f"{payload['wild_type']['forward_fraction']}\t{payload['rnh201']['forward_fraction']}\t"
        f"{strand.get('forward_fraction_difference_rnh201_minus_wt')}\n"
    )
    for base in BASES:
        lines.append(
            f"nucleotide_composition\tbase_{base}\t"
            f"{payload['wild_type']['base_fractions'][base]}\t"
            f"{payload['rnh201']['base_fractions'][base]}\t\n"
        )
    path.write_text("".join(lines))


def write_html(payload: dict, path: Path) -> None:
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en"><head><meta charset="UTF-8">',
        "<title>WT vs rnh201 rNMP comparison</title>",
        "<style>body{font-family:-apple-system,sans-serif;margin:2rem;max-width:960px;}",
        "table{border-collapse:collapse;width:100%;}th,td{border:1px solid #ccc;padding:0.5rem;}",
        ".note{background:#f8f4e8;padding:0.75rem;}</style></head><body>",
        "<h1>Wild type vs RNase H2 (rnh201)</h1>",
        f"<p><em>{_esc(V1_SCOPE)}</em></p>",
    ]
    if payload.get("status") != "compared":
        parts.append(f"<p>Comparison not run. {_esc(payload.get('reason', ''))}</p>")
    else:
        wild = payload["wild_type"]
        mutant = payload["rnh201"]
        parts.append(
            "<table><tr><th>Metric</th><th>Wild type "
            f"({_esc(wild['accession'])})</th><th>rnh201 ({_esc(mutant['accession'])})</th></tr>"
        )
        rows = [
            ("rNMP counts", wild["total_rnmp_5prime_ends"], mutant["total_rnmp_5prime_ends"]),
            ("Unique positions", wild["unique_rnmp_positions"], mutant["unique_rnmp_positions"]),
            ("Calls per Mb", wild["rnmp_calls_per_mb"], mutant["rnmp_calls_per_mb"]),
            ("Forward fraction", wild["forward_fraction"], mutant["forward_fraction"]),
            ("Mean window strand bias", wild["avg_strand_bias"], mutant["avg_strand_bias"]),
        ]
        for label, left, right in rows:
            parts.append(f"<tr><td>{_esc(label)}</td><td>{_esc(left)}</td><td>{_esc(right)}</td></tr>")
        for base in BASES:
            parts.append(
                f"<tr><td>Base fraction {base}</td>"
                f"<td>{_esc(wild['base_fractions'][base])}</td>"
                f"<td>{_esc(mutant['base_fractions'][base])}</td></tr>"
            )
        parts.append("</table>")
        enrich = payload["enrichment"]
        parts.append("<h2>Enrichment (rnh201 / wild type)</h2><ul>")
        parts.append(f"<li>rNMP counts: {_esc(enrich['total_ends_rnh201_over_wt'])}</li>")
        parts.append(f"<li>Unique positions: {_esc(enrich['unique_positions_rnh201_over_wt'])}</li>")
        parts.append(f"<li>Calls per Mb: {_esc(enrich['calls_per_mb_rnh201_over_wt'])}</li>")
        parts.append("</ul>")
        parts.append(f'<p class="note">{_esc(enrich["note"])}</p>')
        parts.append(f'<p class="note">{_esc(payload["mutation_signature_note"])}</p>')
        parts.append(
            "<p>Polymerase assignment: <strong>experimental_not_assigned</strong>. "
            "v1 does not assign replicative polymerases.</p>"
        )
    parts.append(future_work_section())
    parts.append("</body></html>")
    path.write_text("\n".join(parts) + "\n")


def compare_paths(signature_path: Path, hits_path: Path) -> dict:
    signatures = [_load_json(path) for path in _json_files(signature_path)]
    hits = [_load_json(path) for path in _json_files(hits_path)]
    return compare_samples(signatures, _index_hits(hits))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare wild-type and rnh201 rNMP signatures")
    parser.add_argument("--signatures", required=True, help="Signature JSON file or directory")
    parser.add_argument("--hits", required=True, help="Pathway-hits JSON file or directory")
    parser.add_argument("--output", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--html", required=True)
    args = parser.parse_args()
    payload = compare_paths(Path(args.signatures), Path(args.hits))
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n")
    write_table(payload, Path(args.table))
    write_html(payload, Path(args.html))
    print(f"WT vs rnh201 comparison status={payload['status']} → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
