#!/usr/bin/env bash
# Wave 3: gene / excision / overlay analysis only (no bwa, no fetch).
set -euo pipefail
export REPO="${REPO:-/work/app}"
export RESULTS="${RESULTS:-/work/results}"
export REF_DIR="${REF_DIR:-$REPO/data/refs}"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
echo "==> cluster wave3 REPO=$REPO RESULTS=$RESULTS"
python3 "$REPO/scripts/annotate_genes.py"
python3 "$REPO/scripts/cross_study_stats.py" || true
echo "CLUSTER_WAVE3_DONE"
