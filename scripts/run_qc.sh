#!/usr/bin/env bash
# FastQC + fastp.
# trim_5prime is 0 unless the library has a non-genomic prefix (UMI, barcode,
# tagged base) before the genomic sequence. The rNMP coordinate is taken from
# the first genomic base, so that prefix has to come off before alignment.
# The biological 5' end is not quality-trimmed.
set -euo pipefail

threads="${1:?threads}"
qc_dir="${2:?qc_dir}"
fastq_dir="${3:?fastq_dir}"
trim5="${4:-0}"
adapter="${5:-}"
mkdir -p "$qc_dir" trimmed

# trim_front2 is added only when a second read exists. fastp rejects it on
# single-end input, and the v1 ribose-seq runs are single-end.
trim1=()
if [[ "$trim5" != "0" ]]; then
  trim1+=(--trim_front1 "$trim5")
fi
adapter_args=()
if [[ -n "$adapter" ]]; then
  adapter_args+=(--adapter_sequence "$adapter")
fi

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
    --length_required 20 --disable_trim_poly_g \
    "${trim1[@]}" "${adapter_args[@]}"
elif [[ ${#r2[@]} -eq 1 ]]; then
  trim2=()
  if [[ "$trim5" != "0" ]]; then
    trim2+=(--trim_front2 "$trim5")
  fi
  fastp -w "$threads" -i "${r1[0]}" -I "${r2[0]}" \
    -o trimmed/trimmed_R1.fastq.gz -O trimmed/trimmed_R2.fastq.gz \
    --json "$qc_dir/fastp.json" --html "$qc_dir/fastp.html" \
    --length_required 20 --disable_trim_poly_g --detect_adapter_for_pe \
    "${trim1[@]}" "${trim2[@]}" "${adapter_args[@]}"
else
  echo "expected at most one *_R2.fastq.gz, found ${#r2[@]}" >&2
  exit 1
fi
