# hu-rnmp-pipeline

Open-data research prototype for studying **ribonucleotide (rNMP) incorporation**
and **failed DNA repair** in sequencing data from hydroxyurea (HU) / RNR-stress
conditions, with pathway modeling via Reactome and other services.

## What this is

- **Nextflow** pipeline: fetch public SRA/GEO/ENA → QC → align → rNMP-incorporation
  signal modules → pathway enrichment/model drivers
- **Docker** images: pipeline + results browser
- **Harmonia** (`~/harmonia`): bijective multi-omic matrix harmonization leaf
  (held_out — research tooling only, no commercial use/assignment)
- **Reactome / BioModels / PathwayCommons**: pathway models for DNA/RNA repair,
  dNTP biosynthesis, RNR regulation, and ribonucleotide excision repair (RER)
- **Authentik-gated deployment** on `*.flmanbiosci.net` (hwcopeland RKE2 cluster)

## Scope (v1)

- **Public/open data only** (SRA/GEO/ENA accessions). No PHI, no patient BAMs,
  no dbGaP-controlled data in v1.
- **Research prototype** — no clinical claims, no efficacy assertions.
- **Personal repo** (`jonesnoaht`); deployment on flmanbiosci.net is hosting only.

## Repo layout

```
nextflow/           Nextflow pipeline (main.nf, configs, modules)
docker/             Dockerfiles for pipeline + results browser
pathways/           Curated pathway list + model packs + open-data inventory
docs/               Design spec, deploy plan, data-access guide
harmonia/           Harmonia integration (rules, .hm pipelines, matrix bridge)
scripts/            Utility scripts (SRA queries, Reactome fetch, etc.)
manifests/          K8s manifests (theswamp namespace — for iac PR)
data-access/        API key / registration requirements per source
references/         Literature + accession references
.github/workflows/  CI (lint, test, docker build)
```

## Quick start (local Docker)

```bash
# Build pipeline image
docker build -t hu-rnmp-pipeline:dev docker/pipeline/

# Build results browser
docker build -t hu-rnmp-browser:dev docker/browser/

# Run pipeline on a test accession (SRA)
nextflow run main.nf -profile docker --accession SRP123456

# Run results browser
docker run -p 8050:8050 -v $(pwd)/results:/results hu-rnmp-browser:dev
```

## Pathway coverage

See `pathways/pathways.yaml` for the curated list of Reactome + BioModels
pathways relevant to HU / rNMP incorporation / RNA-DNA repair.

## License

Research prototype. Harmonia IP is held_out (NSF/university). No commercial
use or assignment until cleared.
