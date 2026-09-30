#!/usr/bin/env python3
"""Build a small HTML report from pipeline tables."""

from __future__ import annotations

import argparse
import html
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


def _load_json(path: str) -> dict:
    file_path = Path(path)
    if not file_path.exists():
        return {}
    return json.loads(file_path.read_text())


def _esc(value) -> str:
    return html.escape("" if value is None else str(value))


def _metric(label: str, value) -> str:
    return f'<div class="metric">{_esc(label)}: <strong>{_esc(value)}</strong></div>\n'


def build_report(
    accession,
    qc_dir,
    signal,
    mutations,
    mutation_sig,
    enrichment,
    pathway_hits,
    harmonia_matrix,
    harmonia_json,
    output_dir,
) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for src, dst_name in [
        (signal, "rnmp_signal.tsv"),
        (mutations, "rnmp_mutations.tsv"),
        (enrichment, "pathway_enrichment.tsv"),
        (harmonia_matrix, "harmonia_matrix.tsv"),
    ]:
        source = Path(src)
        if source.exists():
            shutil.copy(source, out / dst_name)

    for src, dst_name in [
        (mutation_sig, "mutation_signature.json"),
        (pathway_hits, "pathway_hits.json"),
        (harmonia_json, "harmonia_complement.json"),
    ]:
        source = Path(src)
        if source.exists():
            shutil.copy(source, out / dst_name)

    qc_path = Path(qc_dir)
    if qc_path.exists() and qc_path.is_dir():
        shutil.copytree(qc_path, out / "qc", dirs_exist_ok=True)

    sig = _load_json(mutation_sig)
    hits = _load_json(pathway_hits)
    harm = _load_json(harmonia_json)

    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="UTF-8">',
        f"<title>rNMP report — {_esc(accession)}</title>",
        "<style>",
        "body { font-family: -apple-system, sans-serif; margin: 2rem; max-width: 960px; }",
        "table { border-collapse: collapse; width: 100%; margin: 1rem 0; }",
        "th, td { border: 1px solid #ccc; padding: 0.5rem; text-align: left; }",
        "th { background: #f5f5f5; }",
        ".metric { display: inline-block; margin: 0.5rem; padding: 0.75rem; border: 1px solid #ddd; }",
        ".note { background: #f8f4e8; padding: 0.75rem; }",
        "</style>",
        "</head>",
        "<body>",
        "<h1>rNMP pipeline report</h1>",
        (
            f"<p>Accession: <code>{_esc(accession)}</code> | "
            f"Generated: {_esc(datetime.now(timezone.utc).isoformat())}</p>"
        ),
        "<p><em>Research prototype. Public data only. No clinical claims.</em></p>",
        "<h2>rNMP calls</h2>",
        _metric("Assay", sig.get("assay_type", "")),
        _metric("Technique", sig.get("technique", "")),
        _metric("rNMP counts", sig.get("total_rnmp_5prime_ends", sig.get("total_rnmp_transitions", 0))),
        _metric("Forward", sig.get("forward_ends", "")),
        _metric("Reverse", sig.get("reverse_ends", "")),
        _metric("Unique positions", sig.get("unique_rnmp_positions", "")),
        "<h3>Base counts</h3>",
        "<table><tr><th>Base</th><th>Count</th></tr>",
    ]
    for base, count in sig.get("base_counts", {}).items():
        parts.append(f"<tr><td>{_esc(base)}</td><td>{_esc(count)}</td></tr>")
    parts.append("</table>")

    parts.append('<h3>Per chromosome</h3><table><tr><th>Chrom</th><th>Forward</th><th>Reverse</th></tr>')
    for chrom, counts in sig.get("per_chromosome", {}).items():
        parts.append(
            f"<tr><td>{_esc(chrom)}</td><td>{_esc(counts.get('forward'))}</td>"
            f"<td>{_esc(counts.get('reverse'))}</td></tr>"
        )
    parts.append("</table>")

    parts.append('<p class="note">' + _esc(sig.get("polymerase_note", sig.get("interpretation", ""))) + "</p>")
    if sig.get("experimental"):
        parts.append('<p class="note">WGS signature is experimental. ' + _esc(sig.get("interpretation", "")) + "</p>")

    parts.append("<h2>Signal tests</h2>")
    parts.append(
        "<table><tr><th>Test</th><th>Status</th><th>Statistic</th><th>p</th><th>Note</th></tr>"
    )
    for test in hits.get("signal_tests", []):
        parts.append(
            "<tr>"
            f"<td>{_esc(test.get('id'))}</td>"
            f"<td>{_esc(test.get('status'))}</td>"
            f"<td>{_esc(test.get('statistic'))}</td>"
            f"<td>{_esc(test.get('pvalue'))}</td>"
            f"<td>{_esc(test.get('note'))}</td>"
            "</tr>"
        )
    parts.append("</table>")
    parts.append(
        f"<p>Pathway over-representation: <strong>{_esc(hits.get('pathway_ora', 'unknown'))}</strong>. "
        "Disabled rows are not scores.</p>"
    )

    parts.append("<h2>Harmonia</h2>")
    parts.append(
        f"<p>Status: {_esc(harm.get('status', 'unknown'))}. {_esc(harm.get('reason', ''))}</p>"
    )
    parts.extend(
        [
            "<h2>Files</h2>",
            "<ul>",
            "<li><code>rnmp_signal.tsv</code> — per-position rNMP calls</li>",
            "<li><code>rnmp_mutations.tsv</code> — count summary</li>",
            "<li><code>mutation_signature.json</code></li>",
            "<li><code>pathway_enrichment.tsv</code></li>",
            "<li><code>pathway_hits.json</code></li>",
            "<li><code>harmonia_matrix.tsv</code></li>",
            "<li><code>qc/</code></li>",
            "</ul>",
            "</body></html>",
        ]
    )
    (out / "report.html").write_text("\n".join(parts) + "\n")
    print(f"Report bundle → {out}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the report bundle")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--qc-dir", required=True)
    parser.add_argument("--signal", required=True)
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--mutation-signature", required=True)
    parser.add_argument("--enrichment", required=True)
    parser.add_argument("--pathway-hits", required=True)
    parser.add_argument("--harmonia-matrix", required=True)
    parser.add_argument("--harmonia-json", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    build_report(
        args.accession,
        args.qc_dir,
        args.signal,
        args.mutations,
        args.mutation_signature,
        args.enrichment,
        args.pathway_hits,
        args.harmonia_matrix,
        args.harmonia_json,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
