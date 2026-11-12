# HU-rNMP Pipeline — Design Spec

**Date:** 2026-08-11  
**Author:** Noah T. Jones  
**Repo:** `jonesnoaht/hu-rnmp-pipeline` (personal, public)  
**Deploy:** `hu-rnmp.flmanbiosci.net` (hwcopeland RKE2, Authentik-gated)

**v1 scope (Noah, 2026-11-11).** v1 is wild-type versus RNase H2 (`rnh201`) rNMP mapping on public yeast ribose-seq. It does not assign replicative polymerases. Hydroxyurea is Future work / v2, below.

## 1. Goal

Open-data research prototype that:
1. Compares public wild-type (`RNH201+`) and RNase H2-mutant (`rnh201`) ribose-seq
2. Identifies **ribonucleotide (rNMP) incorporation** positions, strand bias, and base composition
3. Counts **2–5 bp deletions in short tandem repeats** as the mutation signature of rNMPs left in yeast DNA
4. Maps positions to **RNA/DNA repair pathways** (Reactome, BioModels, PathwayCommons) when that test is enabled
5. Uses **Harmonia** for bijective multi-omic matrix harmonization
6. Runs as **Nextflow in Docker** behind **Authentik-gated** results browser on `flmanbiosci.net`

## 2. Scope

- **Public/open data only** (SRA/GEO/ENA). No PHI, no patient BAMs, no dbGaP.
- **Research prototype** — no clinical claims, no efficacy assertions.
- **Personal repo** under `jonesnoaht`. Deployment on flmanbiosci.net is hosting only.
- **Harmonia = held_out IP** (NSF/university). Research tooling only; no commercial use.

## 3. Architecture

```
                    ┌──────────────────────────┐
                    │   GitHub Actions CI       │
                    │   (lint + docker build)    │
                    └──────────┬───────────────┘
                               │ images
                    ┌──────────▼───────────────┐
                    │   Docker Registry         │
                    │   (ghcr.io/jonesnoaht/)   │
                    └──────────┬───────────────┘
                               │
  Internet → Cloudflare → Cilium Gateway → HTTPRoute (Authentik)
                               │
                    ┌──────────▼───────────────┐
                    │  hu-rnmp-browser (Node)    │
                    │  :8050  theswamp           │
                    │  ┌─────────────────────┐  │
                    │  │ Results browser UI  │  │
                    │  │ Run viewer          │  │
                    │  │ Report download     │  │
                    │  └─────────────────────┘  │
                    └──────────┬───────────────┘
                               │ exec
                    ┌──────────▼───────────────┐
                    │  Nextflow pipeline        │
                    │  (pipeline image)          │
                    │                           │
                    │  1. FETCH_SRA (sratoolkit) │
                    │  2. QC_FASTQ (fastqc+fastp)│
                    │  3. ALIGN_BWA (bwa-mem)    │
                    │  4. CALL_RNMP_SIGNAL       │
                    │  5. RNMP_MUTATION_SCAN     │
                    │  6. PATHWAY_ENRICHMENT     │
                    │  7. HARMONIA_JOIN          │
                    │  8. BUILD_REPORT           │
                    └───────────────────────────┘
```

## 4. Pipeline stages

| Stage | Container | Input | Output |
|-------|-----------|-------|--------|
| FETCH_SRA | staphb/sratoolkit | accession | FASTQ pairs |
| QC_FASTQ | staphb/fastqc | FASTQ | trimmed FASTQ + QC reports |
| ALIGN_BWA | staphb/bwa | FASTQ + ref | sorted BAM |
| CALL_RNMP_SIGNAL | hu-rnmp-pipeline:dev | BAM + VCF | rNMP signal TSV + strand bias TSV |
| RNMP_MUTATION_SCAN | hu-rnmp-pipeline:dev | signal + BAM | mutations TSV + signature JSON |
| PATHWAY_ENRICHMENT | hu-rnmp-pipeline:dev | signals + pathways YAML | enrichment TSV + hits JSON |
| HARMONIA_JOIN | hu-rnmp-pipeline:dev | enrichment + rules | harmonized matrix TSV + complement JSON |
| BUILD_REPORT | hu-rnmp-pipeline:dev | all outputs | report_bundle/ (HTML + data) |

## 5. rNMP biology

RNase H2 removes embedded rNMPs by ribonucleotide excision repair (RER). In an
`rnh201` mutant those rNMPs remain in genomic DNA. v1 maps them directly
(ribose-seq) in wild type versus `rnh201`.

The mutation signature of rNMPs left in yeast DNA is a **Top1-dependent deletion
of 2–5 bp inside a short tandem repeat**, not a T>C change (Nick McElhinny et al.
2010, Nat Chem Biol, PMID 20729855; Kim et al. 2011, Science, PMID 21700875).
The pipeline counts those repeat deletions as the headline variant metric.
Raw T>C / A>G counts are a hidden QC field. Homopolymer runs are not counted
as that short-tandem-repeat signature.

v1 does not assign replicative polymerases. Window labels are descriptive
strand-count ratios.

## Future work / v2

Hydroxyurea, with or without RNase H2 and along the Pol η axis, is not a v1
result (Meroni et al. 2019, Nucleic Acids Research, PMID 30847483). The
earlier draft of this spec treated hydroxyurea / RNR-stress sequencing as the
v1 dataset. That contrast is deferred.

Polymerase assignment is a v2 question only for polymerase-variant strains,
relative to replication origins and OK-seq (Daigaku et al. 2015, PMID 25664722).

## 6. Pathway coverage

13 curated pathways in `pathways/pathways.yaml`:

1. RNR regulation (upstream driver)
2. dNTP biosynthesis
3. RNA-DNA damage response / RER (RNase H2)
4. Base excision repair (BER)
5. Mismatch repair (MMR)
6. Nucleotide excision repair (NER)
7. DNA replication stress & fork stalling
8. p53-mediated DNA damage response
9. Homologous recombination & NHEJ
10. Fanconi anemia pathway
11. Cell cycle checkpoints
12. RNA processing / R-loop resolution
13. Oxidative stress response (NRF2/ROS)

Each pathway has Reactome IDs + keywords. Agent team fills model_packs (BioModels SBML,
Reactome FI networks, PathwayCommons SIF).

## 7. Harmonia integration

- `harmonia/rnmp.rules.txt` — Layer-0 rule file for matrix harmonization
- Harmonia compiles from `~/harmonia` (Haskell); binary included in pipeline image
- Used as a **leaf**: harmonize multi-sample enrichment matrices (bijective)
- Records **complements** for lossy transforms (audit trail)
- **held_out IP**: research tooling only; no commercial claims until NSF/university clearance

## 8. Deployment (hwcopeland RKE2)

| Resource | Who applies | Method |
|----------|-------------|--------|
| Deployment + Service + PVC (`theswamp`) | FMB (namespace-admin) | kubectl apply |
| HTTPRoute (`hu-rnmp.flmanbiosci.net`) | Flux/cluster-admin | PR to hwcopeland/iac |
| DNSRecord (A record) | Flux/cluster-admin | PR to iac |
| Authentik provider + outpost | Flux/cluster-admin | PR blueprints + CM |

Pattern mirrors `site-tracker` deployment (`deployment-tracker.yaml`).

## 9. Data access

- **SRA**: NCBI E-utilities; API key recommended (3 req/s → 10 req/s). Free registration.
- **ENA**: Portal API; no key needed.
- **GEO**: FTP/HTTPS bulk download; no key needed.
- **Reactome**: Free downloads (BioPAX, SBML); no key.
- **BioModels**: Free SBML downloads; no key.
- **PathwayCommons**: Free SIF/BIOPAX; no key.
- **KEGG**: REST API free for academic; FTP subscription for bulk.

(Full inventory in `data-access/api_keys_inventory.yaml` — agent-filled.)

## 10. CI/CD

- GitHub Actions: lint Python, build + push Docker images to ghcr.io
- Pipeline image: `ghcr.io/jonesnoaht/hu-rnmp-pipeline:latest`
- Browser image: `ghcr.io/jonesnoaht/hu-rnmp-browser:latest`
- Production: images mirrored to `zot.hwcopeland.net` via Flux image policy

## 11. Next steps

1. Agent team returns with model packs + dataset inventory → review
2. Register NCBI API key (3 req/s is too slow)
3. Build Docker images (local first, then CI)
4. Dry-run Nextflow on the v1 wild-type vs rnh201 ribose-seq profile (`-profile test,docker`)
5. PR to hwcopeland/iac for HTTPRoute + DNS + Authentik
6. Deploy to theswamp
7. First real pipeline run on hu-rnmp.flmanbiosci.net
