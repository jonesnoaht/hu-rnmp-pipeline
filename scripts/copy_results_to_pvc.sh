#!/usr/bin/env bash
# Copy local HU-rNMP results onto the in-cluster PVC via the browser pod.
# PVC is RWO and already mounted at /results on hu-rnmp-browser.
set -euo pipefail
SRC="${1:-/home/noahtjones/hu-rnmp-pipeline/results}"
POD=$(kubectl -n theswamp get pod -l app=hu-rnmp-browser -o jsonpath='{.items[0].metadata.name}')
echo "pod=$POD"
# Copy each top-level run dir (skip lost+found)
for d in "$SRC"/*; do
  [[ -d "$d" ]] || continue
  name=$(basename "$d")
  echo "copy $name → $POD:/results/$name"
  kubectl -n theswamp exec "$POD" -- mkdir -p "/results/$name"
  kubectl -n theswamp cp "$d/." "theswamp/$POD:/results/$name"
done
echo '=== /results after copy ==='
kubectl -n theswamp exec "$POD" -- ls -la /results
