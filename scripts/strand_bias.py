#!/usr/bin/env python3
"""Strand bias: compute replication vs transcription strand coverage.

rNMP incorporation is asymmetric: RNase H2 removes rNMPs, but residual
rNMPs show strand bias. This module computes per-window coverage
asymmetry to identify regions of potential rNMP accumulation.
"""
import argparse
import sys
from collections import defaultdict

import pysam


def compute_strand_bias(bam_path, output_tsv, window=10_000):
    """Compute per-window strand coverage bias."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    results = []
    for chrom in bam.references:
        length = bam.get_reference_length(chrom)
        for start in range(0, length, window):
            end = min(start + window, length)
            fwd = 0
            rev = 0
            for read in bam.fetch(chrom, start, end):
                if read.is_forward:
                    fwd += 1
                else:
                    rev += 1
            total = fwd + rev
            if total < 10:
                continue
            bias = fwd / total
            results.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "fwd_reads": fwd,
                "rev_reads": rev,
                "strand_bias": round(bias, 4),
            })
    bam.close()
    
    with open(output_tsv, "w") as f:
        f.write("chrom\tstart\tend\tfwd_reads\trev_reads\tstrand_bias\n")
        f.writelines(f"{r['chrom']}\t{r['start']}\t{r['end']}\t"
                    f"{r['fwd_reads']}\t{r['rev_reads']}\t{r['strand_bias']}\n" for r in results)
    
    print(f"Strand bias: {len(results)} windows → {output_tsv}", file=sys.stderr)
    return results


def main():
    ap = argparse.ArgumentParser(description="Compute strand bias from BAM")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--window", type=int, default=10_000)
    args = ap.parse_args()
    compute_strand_bias(args.bam, args.output, args.window)


if __name__ == "__main__":
    main()
