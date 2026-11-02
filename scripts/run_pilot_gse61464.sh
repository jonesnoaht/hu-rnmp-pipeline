#!/usr/bin/env bash
# First real HU-rNMP run: Koh 2015 ribose-seq (GSE61464 / PRJNA261234).
# Uses host nix bwa/samtools (home ACLs block docker-as-root on this box).
set -euo pipefail

REPO="${REPO:-/home/noahtjones/hu-rnmp-pipeline}"
PY="$REPO/.venv/bin/python"
REF="$REPO/data/refs/sacCer3.fa"
OUT="$REPO/results/GSE61464"
THREADS="${THREADS:-6}"

mkdir -p "$OUT"

if [[ ! -f "${REF}.bwt" ]]; then
  echo "==> bwa index sacCer3"
  nix-shell -p bwa --run "bwa index '$REF'"
fi

align_one() {
  local acc="$1"
  local dest="$OUT/$acc"
  mkdir -p "$dest"
  if [[ -f "$dest/${acc}.sorted.bam.bai" ]]; then
    echo "==> skip align $acc"
    return 0
  fi
  echo "==> align $acc"
  nix-shell -p bwa samtools --run "
    bwa mem -t $THREADS '$REF' '$REPO/data/raw/GSE61464/${acc}.fastq.gz' \
      | samtools sort -@ $THREADS -o '$dest/${acc}.sorted.bam'
    samtools index '$dest/${acc}.sorted.bam'
    samtools flagstat '$dest/${acc}.sorted.bam' | tee '$dest/flagstat.txt'
  "
}

analyze_one() {
  local acc="$1"
  local dest="$OUT/$acc"
  local bam="$dest/${acc}.sorted.bam"
  echo "==> rNMP signal $acc"
  "$PY" "$REPO/scripts/rnmp_signal.py" \
    --bam "$bam" \
    --output "$dest/rnmp_signal.tsv" \
    --assay ribose_seq \
    --hotspots "$dest/rnmp_hotspots.tsv"
  echo "==> mutation scan $acc"
  "$PY" "$REPO/scripts/rnmp_mutation_scan.py" \
    --signal "$dest/rnmp_signal.tsv" \
    --strand "$dest/rnmp_signal_asymmetry.tsv" \
    --bam "$bam" \
    --output "$dest/rnmp_mutations.tsv" \
    --signature "$dest/mutation_signature.json" \
    --assay ribose_seq
  echo "==> pathway $acc"
  "$PY" "$REPO/scripts/pathway_enrichment.py" \
    --signal "$dest/rnmp_signal.tsv" \
    --mutations "$dest/rnmp_mutations.tsv" \
    --pathways "$REPO/pathways/pathways.yaml" \
    --output "$dest/pathway_enrichment.tsv" \
    --json "$dest/pathway_hits.json"
  mkdir -p "$dest/qc"
  : > "$dest/harmonia_matrix.tsv"
  echo "==> report $acc"
  "$PY" "$REPO/scripts/build_report.py" \
    --accession "$acc" \
    --qc-dir "$dest/qc" \
    --signal "$dest/rnmp_signal.tsv" \
    --mutations "$dest/rnmp_mutations.tsv" \
    --mutation-signature "$dest/mutation_signature.json" \
    --enrichment "$dest/pathway_enrichment.tsv" \
    --pathway-hits "$dest/pathway_hits.json" \
    --harmonia-matrix "$dest/harmonia_matrix.tsv" \
    --output-dir "$dest/report_bundle"
}

for acc in SRR1575897 SRR1575898 SRR1575899 SRR1575900 SRR1575901 SRR1575902 SRR1575903 SRR1575904 SRR1575905; do
  align_one "$acc"
  analyze_one "$acc"
done

cat > "$OUT/README.txt" <<EOF
GSE61464 / PRJNA261234 / SRP047172
Koh et al. 2015 — Global mapping of ribonucleotides embedded in DNA (ribose-seq)
Organism: Saccharomyces cerevisiae
Assay: ribose_seq
Genome: sacCer3
Runs: 9 (all rnh201-background RER-deficient strains)
EOF

echo "PILOT COMPLETE"
find "$OUT" -maxdepth 2 -type f | sort
