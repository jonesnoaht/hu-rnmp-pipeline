#!/usr/bin/env bash
# FastQC + fastp. The 5' end is not cut: direct assays place the rNMP from that end.
set -euo pipefail

threads="${1:?threads}"
qc_dir="${2:?qc_dir}"
fastq_dir="${3:?fastq_dir}"
mkdir -p "$qc_dir" trimmed

shopt -s nullglob
r1=("$fastq_dir"/*_R1.fastq.gz)
r2=("$fastq_dir"/*_R2.fastq.gz)
if [[ ${#r1[@]} -ne 1 ]]; then
  echo "expected one *_R1.fastq.gz in $fastq_dir, found ${#r1[@]}" >&2
  exit 1
fi

fastqc -o "$qc_dir" -t "$threads" "${r1[@]}" "${r2[@]}"

if [[ ${#r2[@]} -eq 0 ]]; then
  fastp -w "$threads" -i "${r1[0]}" -o trimmed/trimmed_R1.fastq.gz \
    --json "$qc_dir/fastp.json" --html "$qc_dir/fastp.html" \
    --length_required 20 --disable_trim_poly_g
elif [[ ${#r2[@]} -eq 1 ]]; then
  fastp -w "$threads" -i "${r1[0]}" -I "${r2[0]}" \
    -o trimmed/trimmed_R1.fastq.gz -O trimmed/trimmed_R2.fastq.gz \
    --json "$qc_dir/fastp.json" --html "$qc_dir/fastp.html" \
    --length_required 20 --disable_trim_poly_g --detect_adapter_for_pe
else
  echo "expected at most one *_R2.fastq.gz, found ${#r2[@]}" >&2
  exit 1
fi
