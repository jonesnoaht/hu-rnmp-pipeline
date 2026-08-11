#!/usr/bin/env python3
"""rNMP Signal: extract ribonucleotide incorporation signatures from BAM/VCF.

rNMP incorporation into genomic DNA produces characteristic signatures:
- T→C transitions on the nascent (replication) strand — embedded rNMPs
  cause misincorporation during replication
- 2-5 bp deletions at rNMP hotspots (poly-N tracts)
- Strand asymmetry: rNMPs on one strand → bias in mutation spectrum

This module extracts these signals from aligned reads + variant calls.
Research only — no clinical claims.
"""
import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict

import pysam


def extract_rnmp_signal(vcf_path, bam_path, output_tsv):
    """Extract per-position rNMP signal from VCF + BAM."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # rNMP signature mutations: T→C (W strand) and A→G (C strand)
    rnmp_transitions = {("T", "C"), ("A", "G")}
    signals = []
    
    try:
        vcf = pysam.VariantFile(vcf_path)
    except Exception:
        # Fallback if no VCF
        vcf = None
    
    if vcf:
        for rec in vcf.fetch():
            ref = rec.ref.upper()
            for alt in rec.alts or []:
                alt = alt.upper()
                if alt == "<NON_REF>" or len(alt) > 1 or len(ref) > 1:
                    continue
                if (ref, alt) in rnmp_transitions:
                    # Check strand from read
                    depth_fwd = 0
                    depth_rev = 0
                    for pileup in bam.pileup(rec.chrom, rec.pos - 1, rec.pos):
                        if pileup.pos == rec.pos - 1:
                            for read in pileup.pileups:
                                if read.is_forward:
                                    depth_fwd += 1
                                else:
                                    depth_rev += 1
                    signals.append({
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": ref,
                        "alt": alt,
                        "type": "rNMP_transition",
                        "depth_fwd": depth_fwd,
                        "depth_rev": depth_rev,
                        "strand_bias": depth_fwd / max(depth_fwd + depth_rev, 1),
                    })
                # Short deletions (2-5 bp) — rNMP slippage
                elif len(ref) >= 2 and len(ref) <= 5 and len(alt) == 1:
                    signals.append({
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": ref,
                        "alt": alt,
                        "type": "rNMP_deletion_slippage",
                        "del_len": len(ref) - 1,
                    })
    bam.close()
    
    # Write TSV
    with open(output_tsv, "w") as f:
        f.write("chrom\tpos\tref\talt\ttype\tdepth_fwd\tdepth_rev\tstrand_bias\tdel_len\n")
        for s in signals:
            f.write("\t".join(str(s.get(k, "")) for k in [
                "chrom", "pos", "ref", "alt", "type",
                "depth_fwd", "depth_rev", "strand_bias", "del_len"
            ]) + "\n")
    
    print(f"rNMP signal: {len(signals)} candidate sites → {output_tsv}", file=sys.stderr)
    return signals


def main():
    ap = argparse.ArgumentParser(description="Extract rNMP signal from VCF+BAM")
    ap.add_argument("--vcf", required=True, help="VCF file (from mpileup/bcftools)")
    ap.add_argument("--bam", required=True, help="Sorted BAM file")
    ap.add_argument("--output", required=True, help="Output TSV")
    args = ap.parse_args()
    extract_rnmp_signal(args.vcf, args.bam, args.output)


if __name__ == "__main__":
    main()
