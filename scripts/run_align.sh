#!/usr/bin/env bash
# End-to-end bowtie2. Local alignment is not used: 5' soft clips move rNMP coordinates.
set -euo pipefail

refdir="${1:?refdir}"
threads="${2:?threads}"
max_insert="${3:?max_insert}"
accession="${4:?accession}"
fastq_dir="${5:?fastq_dir}"

shopt -s nullglob
bt2=("$refdir"/*.1.bt2)
r1=("$fastq_dir"/*_R1.fastq.gz)
r2=("$fastq_dir"/*_R2.fastq.gz)
if [[ ${#bt2[@]} -lt 1 ]]; then
  echo "no bowtie2 index in $refdir" >&2
  exit 1
fi
if [[ ${#r1[@]} -ne 1 ]]; then
  echo "expected one R1 FASTQ, found ${#r1[@]}" >&2
  exit 1
fi
prefix="${bt2[0]%.1.bt2}"

if [[ ${#r2[@]} -eq 0 ]]; then
  bowtie2 --end-to-end --sensitive --seed 1 --no-unal -p "$threads" \
    -x "$prefix" -U "${r1[0]}" 2>align.log \
    | samtools sort -@ "$threads" -o "${accession}.sorted.bam"
elif [[ ${#r2[@]} -eq 1 ]]; then
  bowtie2 --end-to-end --sensitive --seed 1 --no-unal -p "$threads" \
    -x "$prefix" -1 "${r1[0]}" -2 "${r2[0]}" -X "$max_insert" 2>align.log \
    | samtools sort -@ "$threads" -o "${accession}.sorted.bam"
else
  echo "expected at most one R2 FASTQ" >&2
  exit 1
fi
samtools index "${accession}.sorted.bam"
cat align.log
