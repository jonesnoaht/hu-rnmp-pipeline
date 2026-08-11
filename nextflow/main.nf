#!/usr/bin/env nextflow

// HU-rNMP Pipeline: Hydroxyurea + Ribonucleotide Incorporation & DNA Repair
// Open-data prototype — public SRA/GEO/ENA only
// Research use only; no clinical claims.

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────
params.accession      = null   // SRA accession: SRP/GSE/SRR
params.input_csv      = null   // CSV: accession,organism,treatment,control,assay
params.workdir_result = "${params.outdir ?: './results'}"
params.genome         = 'GRCh38'   // also: 'sacCer3' for yeast
params.assay          = 'hyden_seq'  // hyden_seq | ribose_seq | wgs
params.markers_rnmp   = "${params.projectDir}/pathways/rnmp_markers.yaml"
params.pathways       = "${params.projectDir}/pathways/pathways.yaml"
params.harmonia_rules = "${params.projectDir}/harmonia/rnmp.rules.txt"

// ─── Config ───────────────────────────────────────────────────
include { FETCH_SRA;
          QC_FASTQ;
          ALIGN_BWA;
          CALL_RNMP_SIGNAL;
          RNMP_MUTATION_SCAN;
          PATHWAY_ENRICHMENT;
          HARMONIA_JOIN;
          BUILD_REPORT } from './workflows/stages.nf'

// ─── Channel ──────────────────────────────────────────────────
// Input: either --accession (single) or --input_csv (batch)
if (params.accession) {
    accessions_ch = Channel.value([
        [accession: params.accession, organism: params.organism ?: 'human', treatment: 'hydroxyurea', control: params.control ?: 'untreated']
    ])
} else if (params.input_csv) {
    accessions_ch = Channel.fromPath(params.input_csv)
        | splitCsv(header: true)
        | map { row -> [
            accession: row.accession,
            organism: row.organism,
            treatment: row.treatment,
            control: row.control ?: 'untreated'
        ]}
} else {
    log.error "Provide --accession SRP123456 or --input_csv samples.csv"
    System.exit(1)
}

// ─── Pipeline ─────────────────────────────────────────────────

workflow {
    // 1. Fetch raw FASTQ from SRA/ENA
    FETCH_SRA(accessions_ch)

    // 2. QC: FastQC + trimming
    QC_FASTQ(FETCH_SRA.out)

    // 3. Align: BWA-MEM (human) or BWA-MEM (yeast)
    ALIGN_BWA(QC_FASTQ.out, params.genome)

    // 4. rNMP signal: 5' end counting (HydEn-seq) or variant signatures (WGS)
    CALL_RNMP_SIGNAL(ALIGN_BWA.out, params.assay)

    // 5. Mutation scan: analyze rNMP patterns + polymerase strand bias
    RNMP_MUTATION_SCAN(ALIGN_BWA.out, CALL_RNMP_SIGNAL.out, params.assay)

    // 6. Pathway enrichment: map gene hits to Reactome pathway models
    PATHWAY_ENRICHMENT(
        CALL_RNMP_SIGNAL.out,
        RNMP_MUTATION_SCAN.out,
        file(params.pathways)
    )

    // 7. Harmonia: harmonize multi-sample matrices (bijective)
    HARMONIA_JOIN(
        PATHWAY_ENRICHMENT.out,
        file(params.harmonia_rules)
    )

    // 8. Build final report bundle
    BUILD_REPORT(
        QC_FASTQ.out,
        CALL_RNMP_SIGNAL.out,
        RNMP_MUTATION_SCAN.out,
        PATHWAY_ENRICHMENT.out,
        HARMONIA_JOIN.out
    )
}

// ─── Printed summary ──────────────────────────────────────────
workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════╗
    ║  HU-rNMP Pipeline — COMPLETE            ║
    ╠══════════════════════════════════════════╣
    ║  Results: ${params.workdir_result}
    ║  Accession: ${params.accession ?: 'batch CSV'}
    ║  Genome: ${params.genome}
    ╚══════════════════════════════════════════╝
    """.stripIndent()
}
