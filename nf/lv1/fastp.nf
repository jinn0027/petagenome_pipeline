#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastp_fastp_memory = params.memory
params.fastp_fastp_threads = params.threads

params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process fastp {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastp/fastp.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.fastp_fastp_memory}"
    def threads = "${params.fastp_fastp_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        val(p)
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out_{1,2}.fastq", arity: '2')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        fastp \\
            -w ${threads} \\
            --low_complexity_filter \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${pair_id}/out_1.fastq \\
            -O ${pair_id}/out_2.fastq \\
            --cut_front --cut_tail \\
            --cut_mean_quality ${getParam(p, 'fastp_cut_mean_quality')} \\
            --length_required ${getParam(p, 'fastp_reads_minlength')}
        """
}

workflow {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.test_fastp_reads)
    out = fastp(p, reads)
    out.view { i -> "$i" }
}
