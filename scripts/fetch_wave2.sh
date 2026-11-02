#!/usr/bin/env bash
# Wave 2 FASTQs from ENA. SCOPE=wgs (default) | hyden | all
set -euo pipefail
REPO="${REPO:-/home/noahtjones/hu-rnmp-pipeline}"
ROOT="${DATA:-${ROOT:-$REPO/data/raw}}"
SCOPE="${SCOPE:-wgs}"
mkdir -p "$ROOT/GSE125855" "$ROOT/GSE234409"

fetch() {
  local dest="$1" url="$2"
  if [[ -f "$dest" ]] && gzip -t "$dest" 2>/dev/null; then
    echo "ok $(basename "$dest")"
    return 0
  fi
  echo "GET $url"
  wget -c --tries=30 --retry-connrefused --timeout=30 --waitretry=5 -O "$dest" "$url"
  gzip -t "$dest"
}

if [[ "$SCOPE" == "hyden" || "$SCOPE" == "all" ]]; then
  fetch "$ROOT/GSE125855/SRR8503261_1.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/001/SRR8503261/SRR8503261_1.fastq.gz
  fetch "$ROOT/GSE125855/SRR8503261_2.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/001/SRR8503261/SRR8503261_2.fastq.gz
  fetch "$ROOT/GSE125855/SRR8503270_1.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/000/SRR8503270/SRR8503270_1.fastq.gz
  fetch "$ROOT/GSE125855/SRR8503270_2.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/000/SRR8503270/SRR8503270_2.fastq.gz
fi

if [[ "$SCOPE" == "wgs" || "$SCOPE" == "all" ]]; then
  fetch "$ROOT/GSE234409/SRR24871232_1.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR248/032/SRR24871232/SRR24871232_1.fastq.gz
  fetch "$ROOT/GSE234409/SRR24871232_2.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR248/032/SRR24871232/SRR24871232_2.fastq.gz
  fetch "$ROOT/GSE234409/SRR24871236_1.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR248/036/SRR24871236/SRR24871236_1.fastq.gz
  fetch "$ROOT/GSE234409/SRR24871236_2.fastq.gz" \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR248/036/SRR24871236/SRR24871236_2.fastq.gz
fi

echo FETCH_WAVE2_COMPLETE
ls -lh "$ROOT/GSE125855" "$ROOT/GSE234409" || true
