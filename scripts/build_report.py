#!/usr/bin/env python3
"""Build final report bundle: HTML + data files + plots."""
import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


def build_report(accession, qc_dir, signal, mutations, mutation_sig,
                 enrichment, pathway_hits, harmonia_matrix, output_dir):
    """Build self-contained report bundle."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    # Copy data files
    for src, dst_name in [
        (signal, "rnmp_signal.tsv"),
        (mutations, "rnmp_mutations.tsv"),
        (enrichment, "pathway_enrichment.tsv"),
        (harmonia_matrix, "harmonia_matrix.tsv"),
    ]:
        p = Path(src)
        if p.exists():
            shutil.copy(p, out / dst_name)
    
    for src, dst_name in [
        (mutation_sig, "mutation_signature.json"),
        (pathway_hits, "pathway_hits.json"),
    ]:
        p = Path(src)
        if p.exists():
            shutil.copy(p, out / dst_name)
    
    # Copy QC dir if present
    if Path(qc_dir).exists() and Path(qc_dir).is_dir():
        shutil.copytree(qc_dir, out / "qc", dirs_exist_ok=True)
    
    # Load signature for HTML
    sig = {}
    sig_path = Path(mutation_sig)
    if sig_path.exists():
        with open(sig_path) as f:
            sig = json.load(f)
    
    hits = {}
    hits_path = Path(pathway_hits)
    if hits_path.exists():
        with open(hits_path) as f:
            hits = json.load(f)
    
    # Build HTML report
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HU-rNMP Report — {accession}</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 2rem; max-width: 960px; }}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
th, td {{ border: 1px solid #ccc; padding: 0.5rem; text-align: left; }}
th {{ background: #f5f5f5; }}
.metric {{ display: inline-block; margin: 1rem; padding: 1rem;
           border: 1px solid #ddd; border-radius: 4px; }}
.code {{ background: #f8f8f8; padding: 0.5rem; font-family: monospace; }}
</style>
</head>
<body>
<h1>rNMP Pipeline Report</h1>
<p>Accession: <code>{accession}</code> | Generated: {datetime.now(timezone.utc).isoformat()}</p>
<p><em>Research prototype — no clinical claims. Open public data only.</em></p>

<h2>Summary</h2>
<div class="metric">rNMP transitions: <strong>{sig.get('total_rnmp_transitions', 0)}</strong></div>
<div class="metric">rNMP deletions: <strong>{sig.get('total_rnmp_deletions', 0)}</strong></div>
<div class="metric">Strand bias: <strong>{sig.get('avg_strand_bias', 'N/A')}</strong></div>
<div class="metric">Enrichment score: <strong>{sig.get('rnmp_enrichment_score', 'N/A')}</strong></div>

<h3>Transition Breakdown</h3>
<table><tr><th>Type</th><th>Count</th></tr>
"""
    for t, c in sig.get("transition_breakdown", {}).items():
        html += f"<tr><td>{t}</td><td>{c}</td></tr>\n"
    html += "</table>\n"

    html += "<h3>Deletion Size Distribution</h3>\n<table><tr><th>Size (bp)</th><th>Count</th></tr>\n"
    for d, c in sig.get("deletion_size_distribution", {}).items():
        html += f"<tr><td>{d}</td><td>{c}</td></tr>\n"
    html += "</table>\n"

    html += "<h3>Pathway Enrichment</h3>\n<table><tr><th>Pathway</th><th>Score</th><th>Reactome</th></tr>\n"
    for pw in hits.get("pathways", []):
        html += (f"<tr><td>{pw.get('pathway_name', pw.get('pathway_id', ''))}</td>"
                 f"<td>{pw.get('enrichment_score', '')}</td>"
                 f"<td>{', '.join(pw.get('reactome_ids', []))}</td></tr>\n")
    html += """</table>

<h3>Files</h3>
<ul>
<li><code>rnmp_signal.tsv</code> — per-position rNMP signal</li>
<li><code>rnmp_mutations.tsv</code> — mutation summary</li>
<li><code>mutation_signature.json</code> — signature JSON</li>
<li><code>pathway_enrichment.tsv</code> — pathway scores</li>
<li><code>pathway_hits.json</code> — pathway hits</li>
<li><code>harmonia_matrix.tsv</code> — harmonized matrix (Harmonia, held_out)</li>
<li><code>qc/</code> — FastQC reports</li>
</ul>

<footer><p><small>Harmonia IP held_out (NSF/university). Research tooling only — no commercial use.</small></p></footer>
</body></html>"""
    
    (out / "report.html").write_text(html)
    print(f"Report bundle → {out}", file=__import__('sys').stderr)


def main():
    ap = argparse.ArgumentParser(description="Build report bundle")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--mutations", required=True)
    ap.add_argument("--mutation-signature", required=True)
    ap.add_argument("--enrichment", required=True)
    ap.add_argument("--pathway-hits", required=True)
    ap.add_argument("--harmonia-matrix", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    build_report(args.accession, args.qc_dir, args.signal, args.mutations,
                 args.mutation_signature, args.enrichment, args.pathway_hits,
                 args.harmonia_matrix, args.output_dir)


if __name__ == "__main__":
    main()
