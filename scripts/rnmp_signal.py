#!/usr/bin/env python3
"""
rNMP signal extraction from HydEn-seq / ribose-seq / WGS data.

HydEn-seq (Hydrolytic End sequencing, Clausen et al. 2015 NSMB):
  - Alkali cleaves DNA at embedded ribonucleotides (rNMPs)
  - The 5' end of each sequenced read = rNMP incorporation site
  - Strand of the read tells you which strand the rNMP was on
  - In RER-deficient (rnh201Δ) strains, rNMPs accumulate → high signal
  - Polymerase variants (pol1-Y869A, pol2-M644G, pol3-L612G) change
    incorporation patterns → strand-asymmetric rNMP distribution

Ribose-seq (Koh et al. 2015):
  - Similar principle: captures rNMPs via alkaline cleavage
  - Uses rDDP (ribonucleotide-derived DNA product) capture
  - 5' ends map to rNMP positions

For WGS data (indirect):
  - rNMP incorporation leaves mutation signatures (T→C transitions,
    2-5 bp deletions at slippage sites)
  - This is less direct than HydEn-seq/ribose-seq

Key outputs:
  - Per-position rNMP count (5' end pileup)
  - Strand-specific rNMP density
  - Hotspot identification
  - Polymerase-specific incorporation bias

Research only — no clinical claims.
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def extract_hyden_seq_signal(bam_path, output_tsv, window=100):
    """
    Extract rNMP signal from HydEn-seq / ribose-seq BAM.

    The 5' end of each read = rNMP incorporation site.
    Count 5' ends per position, per strand.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Per-position 5' end counts
    forward_ends = defaultdict(int)
    reverse_ends = defaultdict(int)

    # Per-window aggregate (for large genomes)
    window_fwd = defaultdict(int)
    window_rev = defaultdict(int)

    total_reads = 0

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        ref_name = bam.get_reference_name(read.reference_id)

        if read.is_reverse:
            # 5' end of a reverse read = reference position of the last base
            pos = read.reference_end - 1  # 0-based
            reverse_ends[(ref_name, pos)] += 1
            win_start = (pos // window) * window
            window_rev[(ref_name, win_start)] += 1
        else:
            # 5' end of a forward read = reference start
            pos = read.reference_start  # 0-based
            forward_ends[(ref_name, pos)] += 1
            win_start = (pos // window) * window
            window_fwd[(ref_name, win_start)] += 1

        total_reads += 1

    bam.close()

    # Write per-position signal (only non-zero positions)
    with open(output_tsv, "w") as f:
        f.write("chrom\tpos\tstrand\tcount\n")
        for (chrom, pos), count in sorted(forward_ends.items()):
            f.write(f"{chrom}\t{pos}\t+\t{count}\n")
        for (chrom, pos), count in sorted(reverse_ends.items()):
            f.write(f"{chrom}\t{pos}\t-\t{count}\n")

    # Write window summary
    window_tsv = output_tsv.replace(".tsv", "_windows.tsv")
    with open(window_tsv, "w") as f:
        f.write("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tstrand_ratio\n")
        all_keys = set(window_fwd.keys()) | set(window_rev.keys())
        for (chrom, start) in sorted(all_keys):
            fwd = window_fwd.get((chrom, start), 0)
            rev = window_rev.get((chrom, start), 0)
            total = fwd + rev
            ratio = fwd / total if total > 0 else 0.5
            f.write(f"{chrom}\t{start}\t{start + window}\t"
                    f"{fwd}\t{rev}\t{total}\t{ratio:.4f}\n")

    print(f"HydEn-seq signal: {total_reads} reads, "
          f"{len(forward_ends)} fwd 5' ends, "
          f"{len(reverse_ends)} rev 5' ends → {output_tsv}",
          file=sys.stderr)

    return {
        "total_reads": total_reads,
        "forward_ends": len(forward_ends),
        "reverse_ends": len(reverse_ends),
        "window_file": window_tsv,
    }


def extract_wgs_rnmp_signal(vcf_path, bam_path, output_tsv):
    """
    Extract rNMP signal from WGS variant calls (indirect method).

    rNMP incorporation → T→C transitions (replication strand) and
    2-5 bp deletions at poly-N tracts (slippage).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    rnmp_transitions = {("T", "C"), ("A", "G")}
    signals = []

    try:
        vcf = pysam.VariantFile(vcf_path)
    except (OSError, ValueError):
        vcf = None

    if vcf:
        for rec in vcf.fetch():
            ref = rec.ref.upper()
            for alt in rec.alts or []:
                alt = alt.upper()
                if alt == "<NON_REF>" or len(alt) > 1 or len(ref) > 1:
                    continue
                if (ref, alt) in rnmp_transitions:
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
                        "chrom": rec.chrom, "pos": rec.pos, "ref": ref,
                        "alt": alt, "type": "rNMP_transition",
                        "depth_fwd": depth_fwd, "depth_rev": depth_rev,
                        "strand_bias": depth_fwd / max(depth_fwd + depth_rev, 1),
                    })
                elif 2 <= len(ref) <= 5 and len(alt) == 1:
                    signals.append({
                        "chrom": rec.chrom, "pos": rec.pos, "ref": ref,
                        "alt": alt, "type": "rNMP_deletion_slippage",
                        "del_len": len(ref) - 1,
                    })
    bam.close()

    with open(output_tsv, "w") as f:
        f.write("chrom\tpos\tref\talt\ttype\tdepth_fwd\tdepth_rev\tstrand_bias\tdel_len\n")
        f.writelines("\t".join(str(s.get(k, "")) for k in [
                "chrom", "pos", "ref", "alt", "type",
                "depth_fwd", "depth_rev", "strand_bias", "del_len"]) + "\n" for s in signals)

    print(f"WGS rNMP signal: {len(signals)} candidate sites → {output_tsv}",
          file=sys.stderr)
    return {"total_signals": len(signals)}


def compute_strand_asymmetry(bam_path, output_tsv, window=10_000):
    """
    Compute strand-specific rNMP incorporation asymmetry.

    For HydEn-seq: 5' end counts per strand per window.
    For WGS: read coverage per strand per window.

    Strand bias ratio = forward / (forward + reverse)
    A ratio near 0.5 = symmetric (wild-type, RNH201+)
    A ratio near 0.0 or 1.0 = asymmetric (polymerase-specific rNMP incorporation)
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []

    for chrom in bam.references:
        length = bam.get_reference_length(chrom)
        for start in range(0, length, window):
            end = min(start + window, length)
            fwd = 0
            rev = 0
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.is_reverse:
                    rev += 1
                else:
                    fwd += 1
            total = fwd + rev
            if total < 5:
                continue
            bias = fwd / total
            # Classify asymmetry
            if 0.45 <= bias <= 0.55:
                asymmetry = "symmetric"
            elif bias > 0.65:
                asymmetry = "forward_biased"
            elif bias < 0.35:
                asymmetry = "reverse_biased"
            else:
                asymmetry = "mild_bias"
            results.append({
                "chrom": chrom, "start": start, "end": end,
                "fwd_reads": fwd, "rev_reads": rev,
                "strand_bias": round(bias, 4),
                "asymmetry": asymmetry,
            })
    bam.close()

    with open(output_tsv, "w") as f:
        f.write("chrom\tstart\tend\tfwd_reads\trev_reads\tstrand_bias\tasymmetry\n")
        f.writelines(f"{r['chrom']}\t{r['start']}\t{r['end']}\t"
                    f"{r['fwd_reads']}\t{r['rev_reads']}\t{r['strand_bias']}\t"
                    f"{r['asymmetry']}\n" for r in results)

    print(f"Strand asymmetry: {len(results)} windows → {output_tsv}", file=sys.stderr)
    return results


def identify_rnmp_hotspots(signal_tsv, output_tsv, min_count=5, window=500):
    """
    Identify rNMP incorporation hotspots from HydEn-seq 5' end counts.

    Hotspots = windows with significantly higher rNMP counts than
    the genome-wide median. These correspond to regions of high
    polymerase rNMP incorporation or repair failure.
    """
    # Load per-position counts
    counts = defaultdict(lambda: defaultdict(int))
    with open(signal_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, strand, count = parts[0], int(parts[1]), parts[2], int(parts[3])
            counts[chrom][(pos, strand)] += count

    # Aggregate into windows
    window_counts = defaultdict(lambda: {"+": 0, "-": 0})
    for chrom, pos_strands in counts.items():
        for (pos, strand), count in pos_strands.items():
            win_start = (pos // window) * window
            window_counts[(chrom, win_start)][strand] += count

    # Compute median for thresholding
    all_totals = []
    for (chrom, start), strands in window_counts.items():
        total = strands["+"] + strands["-"]
        all_totals.append(total)

    if not all_totals:
        print("No rNMP sites found for hotspot analysis", file=sys.stderr)
        with open(output_tsv, "w") as f:
            f.write("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tfold_over_median\n")
        return []

    all_totals.sort()
    median = all_totals[len(all_totals) // 2]

    # Identify hotspots (fold > 2x median)
    hotspots = []
    for (chrom, start), strands in sorted(window_counts.items()):
        fwd = strands["+"]
        rev = strands["-"]
        total = fwd + rev
        fold = total / median if median > 0 else 0
        if total >= min_count and fold >= 2.0:
            hotspots.append({
                "chrom": chrom, "start": start, "end": start + window,
                "fwd_count": fwd, "rev_count": rev, "total": total,
                "fold_over_median": round(fold, 2),
            })

    with open(output_tsv, "w") as f:
        f.write("chrom\tstart\tend\tfwd_count\trev_count\ttotal\tfold_over_median\n")
        f.writelines(f"{h['chrom']}\t{h['start']}\t{h['end']}\t"
                    f"{h['fwd_count']}\t{h['rev_count']}\t{h['total']}\t"
                    f"{h['fold_over_median']}\n" for h in sorted(hotspots, key=lambda x: x["total"], reverse=True))

    print(f"Hotspots: {len(hotspots)} regions ≥2x median (median={median}) → {output_tsv}",
          file=sys.stderr)
    return hotspots


def main():
    ap = argparse.ArgumentParser(
        description="Extract rNMP signal from HydEn-seq/ribose-seq/WGS data")
    ap.add_argument("--bam", required=True, help="Sorted BAM file")
    ap.add_argument("--output", required=True, help="Output signal TSV")
    ap.add_argument("--assay", default="hyden_seq",
                    choices=["hyden_seq", "ribose_seq", "wgs"],
                    help="Assay type (determines signal extraction method)")
    ap.add_argument("--vcf", default=None,
                    help="VCF file (required for --assay wgs)")
    ap.add_argument("--hotspots", default=None,
                    help="Output hotspots TSV (optional)")
    ap.add_argument("--window", type=int, default=100,
                    help="Window size for aggregation (bp)")
    args = ap.parse_args()

    if args.assay in ("hyden_seq", "ribose_seq"):
        # Direct rNMP detection: 5' ends = rNMP positions
        extract_hyden_seq_signal(args.bam, args.output, args.window)
        if args.hotspots:
            identify_rnmp_hotspots(args.output, args.hotspots, window=500)
    elif args.assay == "wgs":
        # Indirect: variant signatures
        if not args.vcf:
            ap.error("--vcf required for --assay wgs")
        extract_wgs_rnmp_signal(args.vcf, args.bam, args.output)

    # Always compute strand asymmetry
    asymmetry_tsv = args.output.replace(".tsv", "_asymmetry.tsv")
    compute_strand_asymmetry(args.bam, asymmetry_tsv)


if __name__ == "__main__":
    main()
