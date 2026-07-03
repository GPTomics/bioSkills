#!/usr/bin/env nextflow
// Minimal decision-grade RNA-seq pipeline: fastp trim -> Salmon quant -> MultiQC.
// Parses under `nextflow run rnaseq.nf -stub` / `-preview` (wiring only, no tools run).
// Container tags are PINNED to an immutable BioContainers build, never :latest --
// a moving tag breaks reproducibility AND can serve a stale -resume cache hit.
// For production, pin by digest (@sha256:...) so the identity cannot drift.
nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq.gz'
params.salmon_index = 'ref/salmon_index'
params.outdir = 'results'

log.info """
    R N A - S E Q   P I P E L I N E
    ================================
    reads        : ${params.reads}
    salmon_index : ${params.salmon_index}
    outdir       : ${params.outdir}
    """.stripIndent()

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

    stub:
    """
    touch ${sample_id}_trimmed_1.fq.gz ${sample_id}_trimmed_2.fq.gz ${sample_id}_fastp.json
    """
}

process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(index)                                        // shared reference: same index on EVERY sample

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
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir ${sample_id}
    """
}

process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    multiqc . -n multiqc_report
    """

    stub:
    """
    touch multiqc_report.html
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)       // queue: one [id, [r1, r2]] per sample
    index_ch = Channel.fromPath(params.salmon_index, checkIfExists: true)     // queue: ONE item, the shared index

    FASTP(reads_ch)

    // index_ch.first() converts the queue channel to a VALUE channel so the index is reusable
    // on every sample. Passing the bare queue would let sample 1 consume it, leaving the queue
    // empty, and samples 2..N would silently never run (exit 0, no error) -- the #1 dataflow footgun.
    SALMON_QUANT(FASTP.out.reads, index_ch.first())

    qc_ch = FASTP.out.json.collect()                                         // gather all per-sample JSONs
        .mix(SALMON_QUANT.out.quant.map { it[1] }.collect())                 // plus all Salmon output dirs
        .collect()
    MULTIQC(qc_ch)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration: ${workflow.duration}"
    log.info "Success: ${workflow.success}"
}
