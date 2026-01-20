#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "data/*_{1,2}.fq.gz"
params.salmon_index = "ref/salmon_index"
params.outdir = "results"

log.info """
    R N A - S E Q   P I P E L I N E
    ================================
    reads        : ${params.reads}
    salmon_index : ${params.salmon_index}
    outdir       : ${params.outdir}
    """
    .stripIndent()

process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*.fq.gz'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fq.gz"), emit: reads
    path("${sample_id}_fastp.json"), emit: json

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}_trimmed_1.fq.gz \\
        -O ${sample_id}_trimmed_2.fq.gz \\
        --json ${sample_id}_fastp.json \\
        --thread ${task.cpus}
    """
}

process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(index)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant

    script:
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${sample_id} \\
        --threads ${task.cpus} \\
        --validateMappings
    """
}

process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    multiqc . -n multiqc_report
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    index_ch = Channel.fromPath(params.salmon_index, checkIfExists: true)

    FASTP(reads_ch)
    SALMON_QUANT(FASTP.out.reads, index_ch.first())

    qc_ch = FASTP.out.json.collect()
        .mix(SALMON_QUANT.out.quant.map { it[1] }.collect())
        .collect()

    MULTIQC(qc_ch)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration: ${workflow.duration}"
    log.info "Success: ${workflow.success}"
}
