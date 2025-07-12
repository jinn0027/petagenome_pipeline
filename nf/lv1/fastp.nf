#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastp_fastp_memory = params.memory
params.fastp_fastp_threads = params.threads

params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15

process fastp {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastp/fastp.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.fastp_fastp_memory} GB"
    cpus "${params.fastp_fastp_threads}"

    input:
        tuple val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("${pair_id}/out_{1,2}.fastq", arity: '2')
    script:
        """
        mkdir -p ${pair_id}
        fastp \\
            -w ${params.fastp_fastp_threads} \\
            --low_complexity_filter \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${pair_id}/out_1.fastq \\
            -O ${pair_id}/out_2.fastq \\
            --cut_front --cut_tail \\
            --cut_mean_quality ${params.fastp_cut_mean_quality} \\
            --length_required ${params.fastp_reads_minlength}
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_fastp_reads, checkIfExists: true)
    out = fastp(reads)
    out.view { i -> "$i" }
}
