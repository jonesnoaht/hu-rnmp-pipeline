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
# Build pipeline image (context is the repo root; the Dockerfile lives under docker/pipeline)
docker build -f docker/pipeline/Dockerfile -t hu-rnmp-pipeline:dev .

# Build results browser
docker build -t hu-rnmp-browser:dev docker/browser/

# Small public ribose-seq run (SRR1575904, S. cerevisiae, subsampled to 50,000 reads)
nextflow run nextflow/main.nf -profile test,docker

# Run results browser
docker run -p 8050:8050 -v $(pwd)/results:/results hu-rnmp-browser:dev
```

### Reproduce the public-data test

From the repository root, with Docker running and the image tag `hu-rnmp-pipeline:dev`:

```bash
docker build -f docker/pipeline/Dockerfile -t hu-rnmp-pipeline:dev .
nextflow run nextflow/main.nf -profile test,docker
```

`-profile test,docker` reads `samples_test.csv` and keeps 50,000 reads (seed 1) from
[SRR1575904](https://www.ebi.ac.uk/ena/browser/view/SRR1575904) (Koh et al. 2015, PMID 25622106;
S. cerevisiae ribose-seq; single-end MiSeq; 1,068,754 reads; 39,337,368 bytes). The sample is
`pol2-4 rnh201` (KK-107). It is not a hydroxyurea treatment. Reads are aligned end-to-end with
bowtie2 to UCSC sacCer3. Pass `--max_reads 0` on that same command to keep every read.

`samples_hyden_seq.csv` is a corrected study inventory, not the test input. Those series are
not hydroxyurea rNMP maps: SRP182841 and GSE125855 are the same RHII-HydEn-seq study, GSE61464
is single-end ribose-seq, and GSE234409 is a short-treatment WGS copy-number control.

### rNMP coordinate rules

Positions follow Ribose-Map `modules/coordinate.sh` (commit `fff581b`) and Gombolay & Storici,
Nat Protoc 2021, Table 1. Intervals are 0-based and half-open. Paired-end data uses read 1 of
a proper pair only. `hyden_seq` alone is rejected because alkaline and RNase HII HydEn-seq
differ by one nucleotide.

| Technique | Read on + strand | Read on − strand |
|---|---|---|
| ribose-seq | `[start, start+1)` on **−** | `[end-1, end)` on **+** |
| emRiboSeq | `[start-1, start)` on **−** | `[end, end+1)` on **+** |
| Alk-HydEn-seq, Pu-seq | `[start-1, start)` on **+** | `[end, end+1)` on **−** |
| RHII-HydEn-seq | `[start, start+1)` on **+** | `[end-1, end)` on **−** |

Unit tests live in `tests/test_rnmp_coordinates.py` and `tests/test_rnmp_signal.py`.

Pathway over-representation is off unless `--enable_pathway_ora` is set and a gene BED is
supplied. The default stage reports a base-composition chi-square and a strand-balance
binomial test. Harmonia is skipped unless `-params run_harmonia true` and the `harmonia`
binary is on `PATH` (`nextflow run nextflow/main.nf -profile docker --run_harmonia true`).

## Pathway coverage

See `pathways/pathways.yaml` for the curated list of Reactome + BioModels
pathways relevant to HU / rNMP incorporation / RNA-DNA repair.

## License

Research prototype. Harmonia IP is held_out (NSF/university). No commercial
use or assignment until cleared.
