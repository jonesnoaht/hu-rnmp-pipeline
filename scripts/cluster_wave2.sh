#!/usr/bin/env bash
# Cluster entrypoint: fetch + analyze wave-2 on /work (never $HOME).
set -euo pipefail
export REPO="${REPO:-/app}"
export DATA="${DATA:-/work/raw}"
export RESULTS="${RESULTS:-/work/results}"
export PY="${PY:-python3}"
export THREADS="${THREADS:-6}"
export SCOPE="${SCOPE:-wgs}"   # wgs | all
mkdir -p "$DATA/GSE125855" "$DATA/GSE234409" "$RESULTS"
echo "==> cluster wave2 REPO=$REPO DATA=$DATA RESULTS=$RESULTS SCOPE=$SCOPE THREADS=$THREADS"
command -v bwa
command -v samtools
command -v bcftools
command -v python3
REF="${REF:-$REPO/data/refs/sacCer3.fa}"
mkdir -p "$(dirname "$REF")"
if [[ ! -s "$REF" ]]; then
  echo "==> fetch sacCer3"
  wget -O "${REF}.gz" http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  gzip -dc "${REF}.gz" > "$REF"
fi
if [[ ! -s "${REF}.bwt" || ! -s "${REF}.sa" ]]; then
  echo "==> bwa index sacCer3"
  bwa index "$REF"
  samtools faidx "$REF"
fi
ls -lh "$REF"*
"$REPO/scripts/fetch_wave2.sh"
"$REPO/scripts/run_wave2.sh"
"$PY" "$REPO/scripts/cross_study_stats.py" || echo "WARN cross_study_stats failed"
# flatten per-run HTML for the single-segment Express router
for study in GSE125855 GSE234409; do
  shopt -s nullglob
  for html in "$RESULTS/$study"/SRR*/report_bundle/report.html; do
    acc=$(basename "$(dirname "$(dirname "$html")")")
    cp -f "$html" "$RESULTS/$study/${acc}.html"
    echo "flattened $study/${acc}.html"
  done
done
echo CLUSTER_WAVE2_DONE
touch "$RESULTS/CLUSTER_WAVE2_DONE"
ls -lh "$RESULTS" "$RESULTS"/GSE* 2>/dev/null || true
