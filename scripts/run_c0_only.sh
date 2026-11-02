#!/usr/bin/env bash
# Align + call + analyze C0 WGS only (SRR24871236). HU already done.
set -euo pipefail
REPO="${REPO:-/home/noahtjones/hu-rnmp-pipeline}"
# shellcheck disable=SC1091
source /dev/stdin <<FUNCS
$(sed -n '1,96p' "$REPO/scripts/run_wave2.sh")
FUNCS
WRAW="$REPO/data/raw/GSE234409"
WOUT="$REPO/results/GSE234409"
align_pe SRR24871236 "$WRAW/SRR24871236_1.fastq.gz" "$WRAW/SRR24871236_2.fastq.gz" "$WOUT/SRR24871236"
call_wgs SRR24871236 "$WOUT/SRR24871236"
analyze_wgs SRR24871236 "$WOUT/SRR24871236"
echo C0_DONE
