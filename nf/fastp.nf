#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15
params.fastp_memory = 10

process fastp {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastp/fastp.sif"
    publishDir "${params.output}/01_fastp", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}_fastp_out*")
    script:
    """
    fastp -w ${params.threads} -y \\
    -i ${reads[0]} -o ${pair_id}_fastp_out1.fastq.gz \\
    -I ${reads[1]} -O ${pair_id}_fastp_out2.fastq.gz \\
    --cut_front --cut_tail \\
    --cut_mean_quality ${params.fastp_cut_mean_quality} \\
    --length_required ${params.fastp_reads_minlength}
    """
}

workflow {
   raw_short_reads = channel.fromFilePairs(params.reads, checkIfExists: true)
   fastp = fastp(raw_short_reads)
}
