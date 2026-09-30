#!/usr/bin/env nextflow

// Public-data rNMP pipeline. Research use only; no clinical claims.

nextflow.enable.dsl = 2

include { PREPARE_REFERENCE;
          FETCH_ENA;
          QC_FASTQ;
          ALIGN_BOWTIE2;
          CALL_RNMP_SIGNAL;
          RNMP_MUTATION_SCAN;
          PATHWAY_ENRICHMENT;
          HARMONIA_JOIN;
          BUILD_REPORT } from './workflows/stages.nf'

workflow {
    if (!params.accession && !params.input_csv) {
        error "Provide --accession SRR... or --input_csv samples.csv"
    }

    if (params.accession) {
        accessions_ch = Channel.of([
            [
                accession: params.accession,
                organism: params.organism ?: 'Saccharomyces cerevisiae',
                treatment: params.treatment ?: 'unspecified',
                control: params.control ?: 'unspecified',
                assay: params.assay ?: 'ribose_seq',
                technique: params.technique ?: params.assay,
                layout: params.layout ?: ''
            ]
        ])
    } else {
        accessions_ch = Channel.fromPath(params.input_csv, checkIfExists: true)
            | splitCsv(header: true)
            | map { row ->
                [[
                    accession: row.accession,
                    organism: row.organism,
                    treatment: row.treatment,
                    control: row.control ?: 'unspecified',
                    assay: row.assay,
                    technique: row.technique ?: row.assay,
                    layout: row.layout ?: ''
                ]]
            }
    }

    PREPARE_REFERENCE(params.genome)
    ref = PREPARE_REFERENCE.out.dir.first()

    FETCH_ENA(accessions_ch)
    QC_FASTQ(FETCH_ENA.out.reads)
    ALIGN_BOWTIE2(QC_FASTQ.out.clean_fastq, ref)
    CALL_RNMP_SIGNAL(ALIGN_BOWTIE2.out.bam, ref)
    RNMP_MUTATION_SCAN(CALL_RNMP_SIGNAL.out.signal, ref)
    PATHWAY_ENRICHMENT(
        RNMP_MUTATION_SCAN.out.bundle,
        file(params.pathways, checkIfExists: true),
        ref
    )
    HARMONIA_JOIN(
        PATHWAY_ENRICHMENT.out.bundle,
        file(params.harmonia_rules, checkIfExists: true)
    )
    BUILD_REPORT(QC_FASTQ.out.qc.join(HARMONIA_JOIN.out.bundle))
}

workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════╗
    ║  HU-rNMP Pipeline — COMPLETE            ║
    ╠══════════════════════════════════════════╣
    ║  Results: ${params.outdir}
    ║  Genome: ${params.genome}
    ╚══════════════════════════════════════════╝
    """.stripIndent()
}
