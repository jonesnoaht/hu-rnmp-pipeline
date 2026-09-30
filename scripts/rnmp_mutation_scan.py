#!/usr/bin/env python3
"""
rNMP mutation scan: analyze rNMP incorporation patterns.

For HydEn-seq / ribose-seq:
  - 5' end counts per position → rNMP density
  - Strand ratio per window → polymerase identity
  - Hotspot distribution → repair failure regions
  - Compare RNH201+ vs rnh201Δ → RER efficiency

For WGS:
  - T→C transitions (rNMP misincorporation signature)
  - 2-5 bp deletions at poly-N tracts (slippage)
  - Strand-asymmetric mutation burden

Produces:
  - Summary TSV (per-window or per-type counts)
  - Signature JSON (aggregate metrics for pathway enrichment)
"""
import argparse
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path


def scan_hyden_seq(signal_tsv, asymmetry_tsv, output_tsv, signature_json):
    """Analyze HydEn-seq/ribose-seq signal."""
    # Load per-position signal
    pos_counts = []
    strand_counts = {"+": 0, "-": 0}
    with open(signal_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, strand, count = parts[0], int(parts[1]), parts[2], int(parts[3])
            pos_counts.append((chrom, pos, strand, count))
            strand_counts[strand] += count

    # Load asymmetry windows
    asym_windows = []
    with open(asymmetry_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            asym_windows.append({
                "chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                "fwd": int(parts[3]), "rev": int(parts[4]),
                "bias": float(parts[5]), "class": parts[6],
            })

    total_ends = strand_counts["+"] + strand_counts["-"]
    fwd_frac = strand_counts["+"] / max(total_ends, 1)

    # Classify asymmetry distribution
    class_counts = Counter(w["class"] for w in asym_windows)
    bias_values = [w["bias"] for w in asym_windows]
    avg_bias = sum(bias_values) / max(len(bias_values), 1)

    # Per-chromosome rNMP density
    chrom_counts = defaultdict(lambda: {"+": 0, "-": 0})
    for chrom, pos, strand, count in pos_counts:
        chrom_counts[chrom][strand] += count

    # Write summary TSV
    with open(output_tsv, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_5prime_ends\t{total_ends}\n")
        f.write(f"forward_ends\t{strand_counts['+']}\n")
        f.write(f"reverse_ends\t{strand_counts['-']}\n")
        f.write(f"forward_fraction\t{fwd_frac:.4f}\n")
        f.write(f"avg_strand_bias\t{avg_bias:.4f}\n")
        f.write(f"unique_rnmp_positions\t{len(pos_counts)}\n")
        for cls, cnt in class_counts.items():
            f.write(f"window_{cls}\t{cnt}\n")
        f.write("\nchrom\tfwd_ends\trev_ends\ttotal\tratio\n")
        for chrom in sorted(chrom_counts.keys()):
            fwd = chrom_counts[chrom]["+"]
            rev = chrom_counts[chrom]["-"]
            total = fwd + rev
            ratio = fwd / total if total else 0
            f.write(f"{chrom}\t{fwd}\t{rev}\t{total}\t{ratio:.4f}\n")

    # Write signature JSON
    signature = {
        "assay_type": "hyden_seq",
        "total_rnmp_5prime_ends": total_ends,
        "forward_ends": strand_counts["+"],
        "reverse_ends": strand_counts["-"],
        "forward_fraction": round(fwd_frac, 4),
        "unique_rnmp_positions": len(pos_counts),
        "avg_strand_bias": round(avg_bias, 4),
        "strand_asymmetry": ("forward" if fwd_frac > 0.55
                             else "reverse" if fwd_frac < 0.45
                             else "symmetric"),
        "window_classes": dict(class_counts),
        "per_chromosome": {
            chrom: {"forward": chrom_counts[chrom]["+"],
                    "reverse": chrom_counts[chrom]["-"]}
            for chrom in sorted(chrom_counts.keys())
        },
        # RER efficiency proxy: high rNMP count in rnh201Δ vs RNH201+
        # (computed when comparing matched samples; here just raw)
        "rnmp_density_per_mb": round(len(pos_counts) / 12.1, 2),  # yeast ~12.1 Mb
    }

    with open(signature_json, "w") as f:
        json.dump(signature, f, indent=2)

    print(f"HydEn-seq scan: {total_ends} 5' ends, "
          f"{len(pos_counts)} unique rNMP positions → {output_tsv}",
          file=sys.stderr)
    return signature


def scan_wgs(signal_tsv, strand_tsv, output_tsv, signature_json):
    """Analyze WGS rNMP mutation signatures."""
    transitions = []
    deletions = []
    with open(signal_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, mtype = parts[:5]
            if mtype == "rNMP_transition":
                transitions.append({"chrom": chrom, "pos": int(pos),
                                    "ref": ref, "alt": alt})
            elif mtype == "rNMP_deletion_slippage":
                del_len = int(parts[8]) if len(parts) > 8 and parts[8] else 0
                deletions.append({"chrom": chrom, "pos": int(pos),
                                  "del_len": del_len})

    transition_types = Counter(f"{t['ref']}→{t['alt']}" for t in transitions)
    deletion_sizes = Counter(d["del_len"] for d in deletions)
    total = len(transitions) + len(deletions)

    strand_windows = []
    with open(strand_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            strand_windows.append({"chrom": parts[0],
                                   "bias": float(parts[5])})

    bias_values = [w["bias"] for w in strand_windows]
    avg_bias = sum(bias_values) / max(len(bias_values), 1)

    with open(output_tsv, "w") as f:
        f.write("type\tcount\tfrac_of_total\n")
        for mtype, count in transition_types.items():
            f.write(f"transition_{mtype}\t{count}\t{count/max(total,1):.4f}\n")
        for dlen, count in sorted(deletion_sizes.items()):
            f.write(f"deletion_{dlen}bp\t{count}\t{count/max(total,1):.4f}\n")
        f.write(f"TOTAL\t{total}\t1.0\n")

    signature = {
        "assay_type": "wgs",
        "total_rnmp_transitions": len(transitions),
        "total_rnmp_deletions": len(deletions),
        "transition_breakdown": dict(transition_types),
        "deletion_size_distribution": {str(k): v for k, v in deletion_sizes.items()},
        "avg_strand_bias": round(avg_bias, 4),
        "rnmp_enrichment_score": round(len(transitions) / max(total, 1), 4),
    }

    with open(signature_json, "w") as f:
        json.dump(signature, f, indent=2)

    print(f"WGS scan: {len(transitions)} transitions, {len(deletions)} deletions",
          file=sys.stderr)
    return signature


def main():
    ap = argparse.ArgumentParser(description="rNMP mutation scan")
    ap.add_argument("--signal", required=True, help="rNMP signal TSV")
    ap.add_argument("--strand", required=True, help="Strand bias/asymmetry TSV")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--output", required=True, help="Output mutations TSV")
    ap.add_argument("--signature", required=True, help="Output signature JSON")
    ap.add_argument("--assay", default="hyden_seq",
                    choices=["hyden_seq", "ribose_seq", "wgs"])
    args = ap.parse_args()

    if args.assay in ("hyden_seq", "ribose_seq"):
        scan_hyden_seq(args.signal, args.strand, args.output, args.signature)
    else:
        scan_wgs(args.signal, args.strand, args.output, args.signature)


if __name__ == "__main__":
    main()
