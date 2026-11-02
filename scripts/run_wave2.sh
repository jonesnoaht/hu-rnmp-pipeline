#!/usr/bin/env bash
# Wave 2: HydEn-seq pair (GSE125855) + HU/C0 WGS pair (GSE234409).
# Host/Nix OR in-cluster image: uses PATH bwa/samtools/bcftools when present.
set -euo pipefail
REPO="${REPO:-/home/noahtjones/hu-rnmp-pipeline}"
if [[ -n "${PY:-}" ]]; then
  :
elif [[ -x "$REPO/.venv/bin/python" ]]; then
  PY="$REPO/.venv/bin/python"
else
  PY=python3
fi
REF="${REF:-$REPO/data/refs/sacCer3.fa}"
DATA="${DATA:-$REPO/data/raw}"
RESULTS="${RESULTS:-$REPO/results}"
THREADS="${THREADS:-6}"

run_align() {
  local cmd="$1"
  if command -v bwa >/dev/null && command -v samtools >/dev/null; then
    bash -lc "$cmd"
  else
    nix-shell -p bwa samtools --run "$cmd"
  fi
}

run_call() {
  local cmd="$1"
  if command -v bcftools >/dev/null; then
    bash -lc "$cmd"
  else
    nix-shell -p bcftools --run "$cmd"
  fi
}

align_pe() {
  local acc="$1" r1="$2" r2="$3" dest="$4"
  mkdir -p "$dest"
  if [[ -f "$dest/${acc}.sorted.bam.bai" ]]; then
    echo "==> skip align $acc"
    return 0
  fi
  echo "==> align PE $acc"
  run_align "
    bwa mem -t $THREADS '$REF' '$r1' '$r2' \
      | samtools sort -@ $THREADS -o '$dest/${acc}.sorted.bam'
    samtools index '$dest/${acc}.sorted.bam'
    samtools flagstat '$dest/${acc}.sorted.bam' | tee '$dest/flagstat.txt'
  "
}

call_wgs() {
  local acc="$1" dest="$2"
  local bam="$dest/${acc}.sorted.bam"
  local vcf="$dest/${acc}.vcf.gz"
  if [[ -f "$vcf.csi" || -f "$vcf.tbi" ]]; then
    echo "==> skip call $acc"
    return 0
  fi
  echo "==> bcftools call $acc"
  run_call "
    bcftools mpileup -f '$REF' -a FORMAT/AD,FORMAT/DP -Ou '$bam' \
      | bcftools call -mv -Oz -o '$vcf'
    bcftools index '$vcf'
    bcftools stats '$vcf' > '$dest/${acc}.bcftools.stats'
  "
}

analyze_hyden() {
  local acc="$1" dest="$2"
  local bam="$dest/${acc}.sorted.bam"
  echo "==> HydEn signal $acc"
  "$PY" "$REPO/scripts/rnmp_signal.py" \
    --bam "$bam" --output "$dest/rnmp_signal.tsv" \
    --assay hyden_seq --hotspots "$dest/rnmp_hotspots.tsv"
  "$PY" "$REPO/scripts/rnmp_mutation_scan.py" \
    --signal "$dest/rnmp_signal.tsv" \
    --strand "$dest/rnmp_signal_asymmetry.tsv" \
    --bam "$bam" --output "$dest/rnmp_mutations.tsv" \
    --signature "$dest/mutation_signature.json" --assay hyden_seq
  "$PY" "$REPO/scripts/pathway_enrichment.py" \
    --signal "$dest/rnmp_signal.tsv" --mutations "$dest/rnmp_mutations.tsv" \
    --pathways "$REPO/pathways/pathways.yaml" \
    --output "$dest/pathway_enrichment.tsv" --json "$dest/pathway_hits.json"
  mkdir -p "$dest/qc"
  : > "$dest/harmonia_matrix.tsv"
  "$PY" "$REPO/scripts/build_report.py" \
    --accession "$acc" --qc-dir "$dest/qc" \
    --signal "$dest/rnmp_signal.tsv" --mutations "$dest/rnmp_mutations.tsv" \
    --mutation-signature "$dest/mutation_signature.json" \
    --enrichment "$dest/pathway_enrichment.tsv" \
    --pathway-hits "$dest/pathway_hits.json" \
    --harmonia-matrix "$dest/harmonia_matrix.tsv" \
    --output-dir "$dest/report_bundle"
}

analyze_wgs() {
  local acc="$1" dest="$2"
  local bam="$dest/${acc}.sorted.bam"
  local vcf="$dest/${acc}.vcf.gz"
  echo "==> WGS signal $acc"
  "$PY" "$REPO/scripts/rnmp_signal.py" \
    --bam "$bam" --vcf "$vcf" --output "$dest/rnmp_signal.tsv" --assay wgs
  "$PY" "$REPO/scripts/rnmp_mutation_scan.py" \
    --signal "$dest/rnmp_signal.tsv" \
    --strand "$dest/rnmp_signal_asymmetry.tsv" \
    --bam "$bam" --output "$dest/rnmp_mutations.tsv" \
    --signature "$dest/mutation_signature.json" --assay wgs
  "$PY" "$REPO/scripts/pathway_enrichment.py" \
    --signal "$dest/rnmp_signal.tsv" --mutations "$dest/rnmp_mutations.tsv" \
    --pathways "$REPO/pathways/pathways.yaml" \
    --output "$dest/pathway_enrichment.tsv" --json "$dest/pathway_hits.json"
  mkdir -p "$dest/qc"
  : > "$dest/harmonia_matrix.tsv"
  "$PY" "$REPO/scripts/build_report.py" \
    --accession "$acc" --qc-dir "$dest/qc" \
    --signal "$dest/rnmp_signal.tsv" --mutations "$dest/rnmp_mutations.tsv" \
    --mutation-signature "$dest/mutation_signature.json" \
    --enrichment "$dest/pathway_enrichment.tsv" \
    --pathway-hits "$dest/pathway_hits.json" \
    --harmonia-matrix "$dest/harmonia_matrix.tsv" \
    --output-dir "$dest/report_bundle"
}

have_pe() { [[ -f "$1" && -f "$2" ]] && gzip -t "$1" && gzip -t "$2"; }

HRAW="$DATA/GSE125855"
HOUT="$RESULTS/GSE125855"
if have_pe "$HRAW/SRR8503261_1.fastq.gz" "$HRAW/SRR8503261_2.fastq.gz"; then
  align_pe SRR8503261 "$HRAW/SRR8503261_1.fastq.gz" "$HRAW/SRR8503261_2.fastq.gz" "$HOUT/SRR8503261"
  analyze_hyden SRR8503261 "$HOUT/SRR8503261"
else
  echo "==> skip HydEn SRR8503261 (FASTQ missing)"
fi
if have_pe "$HRAW/SRR8503270_1.fastq.gz" "$HRAW/SRR8503270_2.fastq.gz"; then
  align_pe SRR8503270 "$HRAW/SRR8503270_1.fastq.gz" "$HRAW/SRR8503270_2.fastq.gz" "$HOUT/SRR8503270"
  analyze_hyden SRR8503270 "$HOUT/SRR8503270"
else
  echo "==> skip HydEn SRR8503270 (FASTQ missing)"
fi

WRAW="$DATA/GSE234409"
WOUT="$RESULTS/GSE234409"
if have_pe "$WRAW/SRR24871232_1.fastq.gz" "$WRAW/SRR24871232_2.fastq.gz"; then
  align_pe SRR24871232 "$WRAW/SRR24871232_1.fastq.gz" "$WRAW/SRR24871232_2.fastq.gz" "$WOUT/SRR24871232"
  call_wgs SRR24871232 "$WOUT/SRR24871232"
  analyze_wgs SRR24871232 "$WOUT/SRR24871232"
else
  echo "==> skip WGS SRR24871232 (FASTQ missing)"
fi
if have_pe "$WRAW/SRR24871236_1.fastq.gz" "$WRAW/SRR24871236_2.fastq.gz"; then
  align_pe SRR24871236 "$WRAW/SRR24871236_1.fastq.gz" "$WRAW/SRR24871236_2.fastq.gz" "$WOUT/SRR24871236"
  call_wgs SRR24871236 "$WOUT/SRR24871236"
  analyze_wgs SRR24871236 "$WOUT/SRR24871236"
else
  echo "==> skip WGS SRR24871236 (FASTQ missing)"
fi

echo WAVE2_ANALYZE_COMPLETE
