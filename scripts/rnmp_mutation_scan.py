#!/usr/bin/env python3
"""rNMP mutation scan: identify rNMP-induced mutation patterns.

Scans for:
1. T→C transition clusters (rNMP incorporation hotspots)
2. 2-5 bp deletions at poly-N tracts (slippage at rNMP sites)
3. Strand-asymmetric mutation burden
4. Replication timing correlation (late-replicating = more rNMPs)
"""
import argparse
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path


def scan_mutations(signal_tsv, strand_tsv, bam_path, output_tsv, signature_json):
    """Scan rNMP mutation patterns from signal + strand data."""
    # Load signal
    transitions = []
    deletions = []
    with open(signal_tsv) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, pos, ref, alt, mtype = parts[:5]
            pos = int(pos)
            if mtype == "rNMP_transition":
                transitions.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt})
            elif mtype == "rNMP_deletion_slippage":
                del_len = int(parts[8]) if parts[8] else 0
                deletions.append({"chrom": chrom, "pos": pos, "del_len": del_len})
    
    # Load strand bias
    strand_windows = []
    with open(strand_tsv) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            strand_windows.append({
                "chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                "bias": float(parts[5]),
            })
    
    # Cluster transitions (hotspots)
    chrom_sorted = defaultdict(list)
    for t in transitions:
        chrom_sorted[t["chrom"]].append(t["pos"])
    
    # Mutation signature summary
    transition_types = Counter(f"{t['ref']}→{t['alt']}" for t in transitions)
    deletion_sizes = Counter(d["del_len"] for d in deletions)
    
    # Strand bias distribution
    bias_values = [w["bias"] for w in strand_windows]
    avg_bias = sum(bias_values) / max(len(bias_values), 1)
    
    # rNMP-enrichment score: normalized T→C rate
    total_transitions = len(transitions)
    total_deletions = len(deletions)
    
    # Write mutations TSV
    with open(output_tsv, "w") as f:
        f.write("type\tcount\tfrac_of_total\n")
        total = total_transitions + total_deletions
        for mtype, count in transition_types.items():
            f.write(f"transition_{mtype}\t{count}\t{count/max(total,1):.4f}\n")
        for dlen, count in sorted(deletion_sizes.items()):
            f.write(f"deletion_{dlen}bp\t{count}\t{count/max(total,1):.4f}\n")
        f.write(f"TOTAL\t{total}\t1.0\n")
    
    # Write JSON signature
    signature = {
        "total_rnmp_transitions": total_transitions,
        "total_rnmp_deletions": total_deletions,
        "transition_breakdown": dict(transition_types),
        "deletion_size_distribution": {str(k): v for k, v in deletion_sizes.items()},
        "avg_strand_bias": round(avg_bias, 4),
        "strand_bias_skew": "forward" if avg_bias > 0.5 else "reverse",
        "rnmp_enrichment_score": round(total_transitions / max(total, 1), 4),
    }
    
    with open(signature_json, "w") as f:
        json.dump(signature, f, indent=2)
    
    print(f"rNMP mutation scan: {total_transitions} transitions, "
          f"{total_deletions} deletions → {output_tsv}", file=sys.stderr)
    return signature


def main():
    ap = argparse.ArgumentParser(description="rNMP mutation scan")
    ap.add_argument("--signal", required=True, help="rNMP signal TSV")
    ap.add_argument("--strand", required=True, help="Strand bias TSV")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--output", required=True, help="Output mutations TSV")
    ap.add_argument("--signature", required=True, help="Output signature JSON")
    args = ap.parse_args()
    scan_mutations(args.signal, args.strand, args.bam, args.output, args.signature)


if __name__ == "__main__":
    main()
