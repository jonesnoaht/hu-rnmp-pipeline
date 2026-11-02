#!/usr/bin/env bash
# Wait for wave-2 FASTQs (fetch runs separately), analyze HydEn first, then WGS.
set -euo pipefail
REPO=/home/noahtjones/hu-rnmp-pipeline
chmod +x "$REPO/scripts/run_wave2.sh"

wait_files() {
  local timeout="$1"; shift
  local i=0
  while (( i < timeout )); do
    local ok=1
    for f in "$@"; do
      if [[ ! -f "$f" ]] || ! gzip -t "$f" 2>/dev/null; then
        ok=0
        break
      fi
    done
    if (( ok == 1 )); then
      echo "ready: $*"
      return 0
    fi
    sleep 10
    i=$((i + 10))
  done
  echo "timeout waiting for: $*" >&2
  return 1
}

H="$REPO/data/raw/GSE125855"
W="$REPO/data/raw/GSE234409"

echo "==> wait HydEn FASTQs"
wait_files 2400 \
  "$H/SRR8503261_1.fastq.gz" "$H/SRR8503261_2.fastq.gz" \
  "$H/SRR8503270_1.fastq.gz" "$H/SRR8503270_2.fastq.gz"
"$REPO/scripts/run_wave2.sh"
"$REPO/.venv/bin/python" "$REPO/scripts/land_wave2.py" GSE125855

echo "==> wait WGS FASTQs"
wait_files 9000 \
  "$W/SRR24871232_1.fastq.gz" "$W/SRR24871232_2.fastq.gz" \
  "$W/SRR24871236_1.fastq.gz" "$W/SRR24871236_2.fastq.gz"
"$REPO/scripts/run_wave2.sh"
"$REPO/.venv/bin/python" "$REPO/scripts/land_wave2.py" GSE125855 GSE234409

echo ORCH_WAVE2_DONE
