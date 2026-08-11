// Nextflow process modules for HU-rNMP pipeline
// Each process is a stage; called from main.nf

process FETCH_SRA {
    tag "${meta.accession}"
    container 'staphb/sratoolkit:latest'

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("fastq/*_R1.fastq.gz"), path("fastq/*_R2.fastq.gz"), emit: fastq

    script:
    """
    mkdir -p fastq
    prefetch ${meta.accession} -O . 2>/dev/null || true
    fasterqdump ${meta.accession} --outdir fastq --threads ${task.cpus}
    gzip -f fastq/*.fastq 2>/dev/null || true
    # fallback: if paired-end not found, treat as single-end
    ls fastq/ | head -10
    """
}

process QC_FASTQ {
    tag "${meta.accession}"
    container 'staphb/fastqc:latest'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("fastqc/*"), emit: qc
    tuple val(meta), path(r1), path(r2), emit: clean_fastq

    script:
    """
    mkdir -p fastqc
    # FastQC on raw reads
    fastqc -o fastqc -t ${task.cpus} ${r1} ${r2}
    # Fastp for adapter trimming + QC
    fastp -i ${r1} -I ${r2} \
        -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
        --json fastp.json --html fastp.html
    # Use trimmed as output
    cp trimmed_R1.fastq.gz ${r1}
    cp trimmed_R2.fastq.gz ${r2}
    """
}

process ALIGN_BWA {
    tag "${meta.accession}"
    container 'staphb/bwa:latest'

    input:
    tuple val(meta), path(r1), path(r2)
    val genome

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam

    script:
    def ref = genome == 'GRCh38' ? '/refs/GRCh38.fa' : '/refs/sacCer3.fa'
    """
    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o aligned.sorted.bam
    samtools index aligned.sorted.bam
    mv aligned.sorted.bam ${meta.accession}.sorted.bam
    mv aligned.sorted.bam.bai ${meta.accession}.sorted.bam.bai
    """
}

process CALL_RNMP_SIGNAL {
    tag "${meta.accession}"
    container 'hu-rnmp-pipeline:dev'

    input:
    tuple val(meta), path(bam), path(bai)
    val assay

    output:
    tuple val(meta), path("rnmp_signal.tsv"), path("rnmp_signal_windows.tsv"), path("rnmp_signal_asymmetry.tsv"), emit: signal

    script:
    def vcf_flag = assay == "wgs" ? "--vcf raw.vcf" : ""
    def vcf_cmd = assay == "wgs" ? "samtools mpileup -f /refs/genome.fa -vu ${bam} | bcftools call -mv -Ov -o raw.vcf" : ""
    """
    ${vcf_cmd}
    python3 /app/scripts/rnmp_signal.py \
        --bam ${bam} \
        --output rnmp_signal.tsv \
        --assay ${assay} \
        ${vcf_flag} \
        --hotspots rnmp_hotspots.tsv
    """
}

process RNMP_MUTATION_SCAN {
    tag "${meta.accession}"
    container 'hu-rnmp-pipeline:dev'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(signal_tsv), path(windows_tsv), path(asymmetry_tsv)
    val assay

    output:
    tuple val(meta), path("rnmp_mutations.tsv"), path("mutation_signature.json"), emit: mutations

    script:
    """
    python3 /app/scripts/rnmp_mutation_scan.py \
        --signal ${signal_tsv} \
        --strand ${asymmetry_tsv} \
        --bam ${bam} \
        --output rnmp_mutations.tsv \
        --signature mutation_signature.json \
        --assay ${assay}
    """
}

process PATHWAY_ENRICHMENT {
    tag "${meta.accession}"
    container 'hu-rnmp-pipeline:dev'

    input:
    tuple val(meta), path(signal_tsv), path(strand_tsv)
    tuple val(meta2), path(mut_tsv), path(mut_json)
    path pathways_yaml

    output:
    tuple val(meta), path("pathway_enrichment.tsv"), path("pathway_hits.json"), emit: enrichment

    script:
    """
    # Map differentially expressed genes and mutation hotspots
    # to Reactome/BioModels pathways
    python3 /app/scripts/pathway_enrichment.py \
        --signal ${signal_tsv} \
        --mutations ${mut_tsv} \
        --pathways ${pathways_yaml} \
        --model-packs /app/pathways/models/ \
        --output pathway_enrichment.tsv \
        --json pathway_hits.json
    """
}

process HARMONIA_JOIN {
    tag "${meta.accession}"
    container 'hu-rnmp-pipeline:dev'

    input:
    tuple val(meta), path(enrichment_tsv), path(hits_json)
    path rules_file

    output:
    tuple val(meta), path("harmonia_matrix.tsv"), path("harmonia_complement.json"), emit: harmonized

    script:
    """
    # Harmonia: bijective harmonization of multi-sample matrices
    # Uses ~/harmonia Haskell engine (compiled in image)
    # Research tooling only — held_out IP
    harmonia harmonize ${enrichment_tsv} \
        --rules ${rules_file} \
        --output harmonia_matrix.tsv \
        --complement harmonia_complement.json
    """
}

process BUILD_REPORT {
    tag "${meta.accession}"
    container 'hu-rnmp-pipeline:dev'

    input:
    tuple val(meta), path(qc_dir)
    tuple val(meta2), path(signal_tsv), path(strand_tsv)
    tuple val(meta3), path(mut_tsv), path(mut_json)
    tuple val(meta4), path(enrichment_tsv), path(hits_json)
    tuple val(meta5), path(harm_tsv), path(harm_json)

    output:
    path "report_bundle/", emit: report

    script:
    """
    mkdir -p report_bundle
    python3 /app/scripts/build_report.py \
        --accession ${meta.accession} \
        --qc-dir ${qc_dir} \
        --signal ${signal_tsv} \
        --mutations ${mut_tsv} \
        --mutation-signature ${mut_json} \
        --enrichment ${enrichment_tsv} \
        --pathway-hits ${hits_json} \
        --harmonia-matrix ${harm_tsv} \
        --output-dir report_bundle/
    """
}
