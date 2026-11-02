#!/usr/bin/env bash
# Run from the pushed worker image. Uses /app (image) + /work (PVC).
# Does not re-fetch yeast WGS. Does not pull human FASTQ (PVC is 30Gi).
set -euo pipefail

export REPO="${REPO:-/app}"
export RESULTS="${RESULTS:-/work/results}"
export REF_DIR="${REF_DIR:-$REPO/data/refs}"
export DATA="${DATA:-/work/raw}"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export PY="${PY:-python3}"

mkdir -p "$RESULTS/_cross" "$DATA"

echo "==> hu-rnmp image run $(date -Is)"
echo "    REPO=$REPO RESULTS=$RESULTS REF_DIR=$REF_DIR"
echo "    hostname=$(hostname) image_tools:"
command -v bwa
command -v samtools
command -v bcftools
command -v python3
bwa 2>&1 | head -1 || true
samtools --version | head -1
python3 -c "import pysam,scipy,yaml; print('pysam', pysam.__version__, 'scipy', scipy.__version__)"

if [[ ! -s "$REF_DIR/sacCer3.fa.bwt" && -s "$REF_DIR/sacCer3.fa" ]]; then
  echo "==> bwa index sacCer3 (missing .bwt)"
  bwa index "$REF_DIR/sacCer3.fa"
fi

echo "==> gene / excision / FTD-tau overlay"
python3 "$REPO/scripts/annotate_genes.py"

echo "==> cross-study stats (best-effort; missing Koh trees → n=0 is ok)"
python3 "$REPO/scripts/cross_study_stats.py" || true

echo "==> human dataset inventory HTML"
python3 - <<'PY'
import yaml
from pathlib import Path
import os
repo = Path(os.environ.get("REPO", "/app"))
out = Path(os.environ.get("RESULTS", "/work/results")) / "_cross"
src = repo / "catalogs" / "human_datasets.yaml"
data = yaml.safe_load(src.read_text())
recs = data if isinstance(data, list) else data.get("records") or data.get("datasets") or []
if isinstance(data, dict) and not recs:
    for v in data.values():
        if isinstance(v, list) and v and isinstance(v[0], dict) and "accession" in v[0]:
            recs = v
            break
rows = []
for r in recs:
    if not isinstance(r, dict) or not r.get("accession"):
        continue
    rows.append(r)
rows.sort(key=lambda r: ({"CRITICAL":0,"HIGH":1,"MEDIUM":2}.get(str(r.get("relevance","")).upper(), 9), r.get("accession","")))
def cell(x):
    if x is None: return ""
    if isinstance(x, list): return ", ".join(map(str, x))
    return str(x).replace("<","&lt;")
tr = []
for r in rows:
    tr.append("<tr>" + "".join(f"<td>{cell(r.get(k,''))}</td>" for k in
        ("accession","relevance","organism","assay","n_samples","treatment","access")) + "</tr>")
html = f"""<!DOCTYPE html><html><head><meta charset=utf-8>
<title>Human dataset inventory — wave3</title>
<style>
body{{font-family:-apple-system,Segoe UI,sans-serif;margin:1.2rem;max-width:1100px;font-size:14px}}
table{{border-collapse:collapse;width:100%;font-size:12px}}
td,th{{border:1px solid #ccc;padding:.35rem;text-align:left;vertical-align:top}}
th{{background:#f4f4f4}} .warn{{background:#fff6d6;padding:.8rem;border:1px solid #e6d27a}}
</style></head><body>
<h1>Public human / mammalian datasets</h1>
<p>{len(rows)} records from <code>catalogs/human_datasets.yaml</code>. No FASTQ fetched.</p>
<div class="warn">No public human HU WGS that is not sickle-cell. dbGaP = inventory only.
Human RNase H2 disease is AGS, not FTD. Yeast has no MAPT.</div>
<table><thead><tr>
<th>Accession</th><th>Rel</th><th>Organism</th><th>Assay</th><th>n</th><th>Treatment</th><th>Access</th>
</tr></thead><tbody>{''.join(tr)}</tbody></table>
</body></html>"""
out.mkdir(parents=True, exist_ok=True)
(out / "human_datasets.html").write_text(html)
print(f"wrote {out}/human_datasets.html n={len(rows)}")
PY

echo "==> disk"
df -h /work /app | sed -n '1,5p'
echo "CLUSTER_FROM_IMAGE_DONE"
