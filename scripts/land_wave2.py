#!/usr/bin/env python3
"""Refresh study index HTML and copy completed wave-2 runs onto the PVC."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

REPO = Path("/home/noahtjones/hu-rnmp-pipeline")
STUDIES = {
    "GSE125855": {
        "title": "GSE125855 — HydEn-seq (Zhou / Kunkel)",
        "blurb": "S. cerevisiae rnh201 HydEn-seq: NT vs RNase HII. Genome sacCer3. Research prototype — no clinical claims.",
        "runs": [
            ("SRR8503261", "rnh201 NT (TAK_1811)"),
            ("SRR8503270", "rnh201 RNaseHII (TAK_1827)"),
        ],
    },
    "GSE234409": {
        "title": "GSE234409 — HU vs C0 WGS pair",
        "blurb": "Matched 2-1 series: HU-14-B vs untreated C0-6. Yeast WGS rNMP-proxy variants. Genome sacCer3. Research prototype — no clinical claims.",
        "runs": [
            ("SRR24871232", "2-1-HU-14-B"),
            ("SRR24871236", "2-1-C0-6"),
        ],
    },
}


def pod_name() -> str:
    return subprocess.check_output(
        [
            "kubectl",
            "-n",
            "theswamp",
            "get",
            "pod",
            "-l",
            "app=hu-rnmp-browser",
            "-o",
            "jsonpath={.items[0].metadata.name}",
        ],
        text=True,
    ).strip()


def write_index(study: str, spec: dict) -> Path:
    root = REPO / "results" / study
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for acc, label in spec["runs"]:
        sigp = root / acc / "mutation_signature.json"
        htmlp = root / f"{acc}.html"
        src = root / acc / "report_bundle" / "report.html"
        if src.exists():
            shutil.copy2(src, htmlp)
        if sigp.exists() and htmlp.exists():
            sig = json.loads(sigp.read_text())
            rows.append(
                (
                    acc,
                    label,
                    sig.get("total_rnmp_5prime_ends", sig.get("total_rnmp_transitions", "—")),
                    sig.get("unique_rnmp_positions", sig.get("total_rnmp_deletions", "—")),
                    sig.get("rnmp_density_per_mb", sig.get("rnmp_enrichment_score", "—")),
                    sig.get("strand_asymmetry", sig.get("avg_strand_bias", "—")),
                    True,
                )
            )
        else:
            rows.append((acc, label, "…", "…", "…", "running", False))
    html = [
        "<!DOCTYPE html><html><head><meta charset=utf-8>",
        f"<title>{spec['title']}</title>",
        "<style>body{font-family:-apple-system,sans-serif;margin:2rem;max-width:980px}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccc;padding:.45rem;text-align:left}th{background:#f5f5f5}</style></head><body>",
        f"<h1>{spec['title']}</h1>",
        f"<p>{spec['blurb']}</p>",
        "<table><tr><th>Run</th><th>Sample</th><th>Signal</th><th>Sites / dels</th><th>Density / score</th><th>Strand</th></tr>",
    ]
    for acc, label, a, b, c, d, ok in rows:
        cell = f'<a href="{acc}.html">{acc}</a>' if ok else acc
        html.append(f"<tr><td>{cell}</td><td>{label}</td><td>{a}</td><td>{b}</td><td>{c}</td><td>{d}</td></tr>")
    html.append("</table></body></html>")
    (root / "report.html").write_text("\n".join(html))
    (root / "index.html").write_text("\n".join(html))
    return root


def land(study: str) -> None:
    spec = STUDIES[study]
    root = write_index(study, spec)
    pod = pod_name()
    subprocess.check_call(["kubectl", "-n", "theswamp", "exec", pod, "--", "mkdir", "-p", f"/results/{study}"])
    for acc, _ in spec["runs"]:
        src = root / acc
        if (src / "report_bundle" / "report.html").exists():
            subprocess.check_call(
                ["kubectl", "-n", "theswamp", "exec", pod, "--", "mkdir", "-p", f"/results/{study}/{acc}"]
            )
            subprocess.check_call(
                [
                    "kubectl",
                    "-n",
                    "theswamp",
                    "cp",
                    f"{src}/.",
                    f"{pod}:/results/{study}/{acc}",
                    "--no-preserve=true",
                ]
            )
            html = root / f"{acc}.html"
            if html.exists():
                subprocess.check_call(
                    [
                        "kubectl",
                        "-n",
                        "theswamp",
                        "cp",
                        str(html),
                        f"{pod}:/results/{study}/{acc}.html",
                        "--no-preserve=true",
                    ]
                )
    subprocess.check_call(
        [
            "kubectl",
            "-n",
            "theswamp",
            "cp",
            str(root / "report.html"),
            f"{pod}:/results/{study}/report.html",
            "--no-preserve=true",
        ]
    )
    print(f"landed {study}")


if __name__ == "__main__":
    import sys

    targets = sys.argv[1:] or list(STUDIES)
    for s in targets:
        land(s)
