// Nextflow stages for the public-data rNMP pipeline.

process PREPARE_REFERENCE {
    tag "${genome}"

    input:
    val genome

    output:
    path "reference", emit: dir

    script:
    """
    bash ${params.tool_scripts}/prepare_reference.sh ${genome} reference ${task.cpus}
    """
}

process FETCH_ENA {
    tag "${meta.accession}"

    input:
    val meta

    output:
    tuple val(meta), path("fastq"), emit: reads
    path "fetch_manifest.json", emit: manifest

    script:
    """
    mkdir -p fastq
    python3 ${params.tool_scripts}/fetch_ena.py \
        --accession ${meta.accession} \
        --outdir fastq \
        --max-reads ${params.max_reads} \
        --seed ${params.subsample_seed} \
        --max-runs ${params.max_runs} \
        --manifest fetch_manifest.json
    """
}

process QC_FASTQ {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("qc"), emit: qc
    tuple val(meta), path("trimmed"), emit: clean_fastq

    script:
    """
    bash ${params.tool_scripts}/run_qc.sh ${task.cpus} qc ${fastq_dir}
    """
}

process ALIGN_BOWTIE2 {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(fastq_dir)
    path refdir

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam

    script:
    """
    bash ${params.tool_scripts}/run_align.sh ${refdir} ${task.cpus} ${params.max_insert} ${meta.accession} ${fastq_dir}
    """
}

process CALL_RNMP_SIGNAL {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(bam), path(bai)
    path refdir

    output:
    tuple val(meta), path("rnmp_signal.tsv"), path("rnmp_signal_windows.tsv"), path("rnmp_signal_asymmetry.tsv"), path("rnmp_signal_stats.json"), emit: signal
    path "rnmp_hotspots.tsv", emit: hotspots

    script:
    def fasta = "\$(ls ${refdir}/*.fa | head -n 1)"
    if (meta.assay == "wgs")
        """
        fasta=${fasta}
        bcftools mpileup -f \$fasta -Ou ${bam} | bcftools call -mv -Ov -o raw.vcf
        python3 ${params.tool_scripts}/rnmp_signal.py \
            --bam ${bam} \
            --fasta \$fasta \
            --assay wgs \
            --vcf raw.vcf \
            --output rnmp_signal.tsv \
            --hotspots rnmp_hotspots.tsv \
            --min-mapq ${params.min_mapq}
        touch rnmp_hotspots.tsv
        """
    else
        """
        fasta=${fasta}
        python3 ${params.tool_scripts}/rnmp_signal.py \
            --bam ${bam} \
            --fasta \$fasta \
            --assay ${meta.assay} \
            --technique '${meta.technique}' \
            --output rnmp_signal.tsv \
            --hotspots rnmp_hotspots.tsv \
            --min-mapq ${params.min_mapq}
        touch rnmp_hotspots.tsv
        """
}

process RNMP_MUTATION_SCAN {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json)
    path refdir

    output:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path("rnmp_mutations.tsv"), path("mutation_signature.json"), emit: bundle

    script:
    def fai = "\$(ls ${refdir}/*.fa.fai | head -n 1)"
    """
    fai=${fai}
    python3 ${params.tool_scripts}/rnmp_mutation_scan.py \
        --signal ${signal_tsv} \
        --strand ${asymmetry_tsv} \
        --stats ${stats_json} \
        --fai \$fai \
        --output rnmp_mutations.tsv \
        --signature mutation_signature.json \
        --assay ${meta.assay} \
        --technique '${meta.technique}'
    """
}

process PATHWAY_ENRICHMENT {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path(mut_tsv), path(sig_json)
    path pathways_yaml
    path refdir

    output:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path(mut_tsv), path(sig_json), path("pathway_enrichment.tsv"), path("pathway_hits.json"), emit: bundle

    script:
    def fasta = "\$(ls ${refdir}/*.fa | head -n 1)"
    def ora = params.enable_pathway_ora ? "--enable-pathway-ora" : ""
    """
    fasta=${fasta}
    python3 ${params.tool_scripts}/pathway_enrichment.py \
        --signal ${signal_tsv} \
        --mutations ${mut_tsv} \
        --pathways ${pathways_yaml} \
        --fasta \$fasta \
        ${ora} \
        --output pathway_enrichment.tsv \
        --json pathway_hits.json
    """
}

process HARMONIA_JOIN {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path(mut_tsv), path(sig_json), path(enrichment_tsv), path(hits_json)
    path rules_file

    output:
    tuple val(meta), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path(mut_tsv), path(sig_json), path(enrichment_tsv), path(hits_json), path("harmonia_matrix.tsv"), path("harmonia_complement.json"), emit: bundle

    script:
    """
    python3 ${params.tool_scripts}/harmonia_join.py \
        --enrichment ${enrichment_tsv} \
        --rules ${rules_file} \
        --output harmonia_matrix.tsv \
        --complement harmonia_complement.json \
        --run ${params.run_harmonia}
    """
}

process BUILD_REPORT {
    tag "${meta.accession}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(qc_dir), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv), path(stats_json), path(mut_tsv), path(sig_json), path(enrichment_tsv), path(hits_json), path(harm_tsv), path(harm_json)

    output:
    path "report_bundle/", emit: report

    script:
    """
    python3 ${params.tool_scripts}/build_report.py \
        --accession ${meta.accession} \
        --qc-dir ${qc_dir} \
        --signal ${signal_tsv} \
        --mutations ${mut_tsv} \
        --mutation-signature ${sig_json} \
        --enrichment ${enrichment_tsv} \
        --pathway-hits ${hits_json} \
        --harmonia-matrix ${harm_tsv} \
        --harmonia-json ${harm_json} \
        --output-dir report_bundle/
    """
}
