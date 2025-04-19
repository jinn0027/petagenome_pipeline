#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15

process fastp {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastp/fastp.sif"
    publishDir "${params.output}/fastp/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}_fastp_out*")
    script:
    """
    fastp \\
        -w ${params.threads} \\
        --low_complexity_filter \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${pair_id}_fastp_out1.fastq.gz \\
        -O ${pair_id}_fastp_out2.fastq.gz \\
        --cut_front --cut_tail \\
        --cut_mean_quality ${params.fastp_cut_mean_quality} \\
        --length_required ${params.fastp_reads_minlength}
    """
}

workflow {
    reads = channel.fromFilePairs(params.test_fastp_reads, checkIfExists: true)
    fastp = fastp(reads)
    //fastp.view { i -> "$i" }
}
