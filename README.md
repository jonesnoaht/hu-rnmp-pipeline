# hu-rnmp-pipeline

Open-data research prototype for **ribonucleotide (rNMP) maps** in yeast.
v1 compares public wild-type (`RNH201+`) ribose-seq with an RNase H2 mutant
(`rnh201`). Pathway models come from Reactome and related public resources.

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

- **Wild type vs RNase H2.** The test profile maps rNMPs in one public
  `RNH201+` library and one matched `rnh201` library (same study, strain, and
  restriction-enzyme set).
- **Public/open data only** (SRA/GEO/ENA accessions). No PHI, no patient BAMs,
  no dbGaP-controlled data in v1. No credentials beyond an optional `NCBI_API_KEY`.
- **Research prototype** — no clinical claims, no efficacy assertions.
- **No polymerase assignment.** `polymerase_assignment` stays
  `experimental_not_assigned`.
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

# Wild type vs rnh201 ribose-seq (two public runs, each subsampled to 50,000 reads)
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

`-profile test,docker` reads `samples_test.csv` and keeps 50,000 reads (seed 1) from each run.
Both libraries are nuclear ribose-seq from Balachander et al. 2020
([PMID 32415081](https://pubmed.ncbi.nlm.nih.gov/32415081/),
[PRJNA613920](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA613920),
[ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA613920)):

| Role | Run | Library | Strain | Reads | FASTQ bytes |
|---|---|---|---|---:|---:|
| Wild type (`RNH201+`) | [SRR11364941](https://www.ebi.ac.uk/ena/browser/view/SRR11364941) | E134 WT RE3 FS162 | E134 | 5,278,518 | 335,187,453 |
| `rnh201` | [SRR11364929](https://www.ebi.ac.uk/ena/browser/view/SRR11364929) | E134 rnh201 RE3 FS116 | E134 | 5,622,971 | 157,309,501 |

Koh et al. 2015 ([PMID 25622106](https://pubmed.ncbi.nlm.nih.gov/25622106/),
[PRJNA261234](https://www.ebi.ac.uk/ena/browser/view/PRJNA261234) /
[GSE61464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61464)) has no `RNH201+` library.
[SRR1575904](https://www.ebi.ac.uk/ena/browser/view/SRR1575904) is `pol2-4 rnh201`, so it is not the v1 contrast.
The 2020 reads start with 12 non-genomic bases (UMI, barcode, and tagged nucleotide) and are
trimmed before alignment, along with the 3' adapter `AGTTGCGACACGGATCTATCA`. UMIs are not collapsed.
Reads are aligned end-to-end with bowtie2 to UCSC sacCer3. The comparison report is
`wt_vs_rnh201/comparison_report.html` (position counts, strand bias, nucleotide composition,
and rnh201/wild-type enrichment). Pass `--max_reads 0` on that same command to keep every read.

`samples_hyden_seq.csv` is a corrected study inventory, not the test input. SRP182841 and
GSE125855 are the same RHII-HydEn-seq study, GSE61464 is the Koh 2015 ribose-seq series
(no wild-type run), and GSE234409 is a WGS copy-number control held for Future work / v2.

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

Unit tests live in `tests/test_rnmp_coordinates.py`, `tests/test_rnmp_signal.py`,
`tests/test_repeat_deletion.py`, and `tests/test_rnmp_compare.py`.

### Mutation signature

For variant calls, the headline metric is a **deletion of 2–5 bp that sits in a short
tandem repeat** (repeat unit 2–5 bp, at least two copies; homopolymer units are not
counted). That is the Top1-dependent signature of rNMPs left in yeast DNA
(Nick McElhinny et al. 2010, PMID 20729855; Kim et al. 2011, PMID 21700875).
Raw T>C and A>G counts stay in `mutation_signature.json` under `qc` and are not
shown as an rNMP signal. Direct ribose-seq runs report rNMP positions instead of
that deletion count.

Pathway over-representation is off unless `--enable_pathway_ora` is set and a gene BED is
supplied. The default stage reports a base-composition chi-square and a strand-balance
binomial test. Harmonia is skipped unless `-params run_harmonia true` and the `harmonia`
binary is on `PATH` (`nextflow run nextflow/main.nf -profile docker --run_harmonia true`).

## Pathway coverage

See `pathways/pathways.yaml` for the curated Reactome + BioModels list used with
rNMP maps. Hydroxyurea wording in that file is background for v2, not a v1 result.

## Future work / v2

Hydroxyurea, with or without RNase H2 and along the Pol η axis, is not a v1 result
(Meroni et al. 2019, Nucleic Acids Research, PMID 30847483).

Polymerase assignment is not done in v1. Assignment is a v2 question only for
polymerase-variant strains, relative to replication origins and OK-seq
(Daigaku et al. 2015, PMID 25664722).

## License

Research prototype. Harmonia IP is held_out (NSF/university). No commercial
use or assignment until cleared.
