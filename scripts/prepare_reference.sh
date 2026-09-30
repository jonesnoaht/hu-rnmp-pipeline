#!/usr/bin/env bash
# Download a public yeast reference and build a bowtie2 index.
set -euo pipefail

genome="${1:?genome}"
outdir="${2:?outdir}"
threads="${3:-2}"
mkdir -p "$outdir"

case "$genome" in
  sacCer3)
    url="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz"
    curl -fL --retry 3 --retry-delay 2 -o "$outdir/sacCer3.fa.gz" "$url"
    gzip -dc "$outdir/sacCer3.fa.gz" > "$outdir/sacCer3.fa"
    rm -f "$outdir/sacCer3.fa.gz"
    bowtie2-build --threads "$threads" "$outdir/sacCer3.fa" "$outdir/sacCer3"
    samtools faidx "$outdir/sacCer3.fa"
    ;;
  *)
    echo "No reference bundle is configured for genome '$genome'. The test profile uses sacCer3." >&2
    exit 1
    ;;
esac
