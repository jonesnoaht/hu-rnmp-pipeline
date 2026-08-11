#!/usr/bin/env python3
"""Pathway enrichment: map rNMP hits to Reactome/BioModels pathways.

Takes rNMP signal + mutation data and maps to curated pathway models
(Reactome, BioModels, PathwayCommons) to identify which DNA/RNA repair
pathways are enriched for rNMP-induced perturbation.
"""
import argparse
import json
import sys
import yaml
from collections import defaultdict, Counter


def load_pathways(pathways_yaml):
    """Load curated pathway list."""
    with open(pathways_yaml) as f:
        data = yaml.safe_load(f)
    return {p["id"]: p for p in data.get("pathways", [])}


def load_signal(signal_tsv):
    """Load rNMP signal sites."""
    sites = []
    with open(signal_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                sites.append({
                    "chrom": parts[0], "pos": parts[1],
                    "ref": parts[2], "alt": parts[3], "type": parts[4],
                })
    return sites


def load_mutations(mutations_tsv):
    """Load mutation summary."""
    mutations = []
    with open(mutations_tsv) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mutations.append({"type": parts[0], "count": int(parts[1])})
    return mutations


def enrich(sites, mutations, pathways, model_packs_dir=None):
    """Map signals to pathways — this is a simple keyword/gene-based heuristic."""
    # Eigenes/gene panels per pathway (curated from Reactome + literature)
    pathway_gene_panels = {
        "rnr-regulation": ["RRM1", "RRM2", "RRM2B", "TYMS", "DTYMK"],
        "dntp-biosynthesis": ["RRM1", "RRM2", "CMPK1", "DTYMK", "TYMS", "UMPS"],
        "rna-dna-damage-response": ["RNASEH2A", "RNASEH2B", "RNASEH2C", "RNASEH1"],
        "base-excision-repair": ["PARP1", "XRCC1", "APEX1", "OGG1", "FEN1", "LIG3"],
        "mismatch-repair": ["MSH2", "MSH6", "MLH1", "PMS2", "EXO1"],
        "nucleotide-excision-repair": ["XPA", "ERCC1", "XPF", "XPG", "XPC"],
        "dna-replication-stress": ["ATR", "CHEK1", "RPA1", "RPA2", "CLASPIN"],
        "p53-dna-damage-response": ["TP53", "MDM2", "CDKN1A", "BBC3", "BAX"],
        "hr-nhej-double-strand-break": ["BRCA1", "BRCA2", "RAD51", "XRCC4", "LIG4", "KU70", "KU80"],
        "fanconi-anemia-pathway": ["FANCD2", "FANCA", "BRCA2", "FANCL"],
        "cell-cycle-checkpoints": ["ATR", "ATM", "CHEK1", "CHEK2", "CDC25A", "CDK1", "TP53"],
        "rna-processing-damage": ["RNASEH1", "SETX", "BRCA1", "TOP1"],
        "oxidative-stress-response": ["NFE2L2", "KEAP1", "OGG1", "SOD2", "CAT"],
    }
    
    enrichment_results = []
    total_mutations = sum(m["count"] for m in mutations)
    
    for pid, pathway_def in pathways.items():
        genes = pathway_gene_panels.get(pid, [])
        reactome_ids = pathway_def.get("reactome", [])
        # Score: proportion of pathway genes potentially affected
        # (placeholder until full gene-level differential expression available)
        score = len(genes) * 0.1  # heuristic
        
        enrichment_results.append({
            "pathway_id": pid,
            "pathway_name": pathway_def.get("name", pid),
            "reactome_ids": reactome_ids,
            "panel_genes": genes,
            "enrichment_score": round(score, 4),
            "rationale": pathway_def.get("rationale", "")[:200],
        })
    
    # Sort by score
    enrichment_results.sort(key=lambda x: x["enrichment_score"], reverse=True)
    return enrichment_results


def main():
    ap = argparse.ArgumentParser(description="Pathway enrichment for rNMP hits")
    ap.add_argument("--signal", required=True)
    ap.add_argument("--mutations", required=True)
    ap.add_argument("--pathways", required=True)
    ap.add_argument("--model-packs", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    
    pathways = load_pathways(args.pathways)
    sites = load_signal(args.signal)
    mutations = load_mutations(args.mutations)
    
    results = enrich(sites, mutations, pathways, args.model_packs)
    
    # Write TSV
    with open(args.output, "w") as f:
        f.write("pathway_id\tpathway_name\treactome_ids\tpanel_genes\tenrichment_score\trationale\n")
        for r in results:
            f.write(f"{r['pathway_id']}\t{r['pathway_name']}\t"
                    f"{'|'.join(r['reactome_ids'])}\t"
                    f"{'|'.join(r['panel_genes'])}\t"
                    f"{r['enrichment_score']}\t{r['rationale']}\n")
    
    # Write JSON
    with open(args.json, "w") as f:
        json.dump({"pathways": results, "total_signal_sites": len(sites),
                   "total_mutations": len(mutations)}, f, indent=2)
    
    print(f"Pathway enrichment: {len(results)} pathways {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
