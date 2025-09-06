#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqc_fastqc_memory = params.memory
params.fastqc_fastqc_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process fastqc {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastqc/fastqc.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.fastqc_fastqc_memory}"
    def threads = "${params.fastqc_fastqc_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        echo "${processProfile(task)}"
        export XDG_CACHE_HOME=\$(pwd)/.cache
        mkdir -p .cache
        mkdir -p ${pair_id}
        fastqc \\
            -t ${threads} \\
            -o ${pair_id} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    reads = createPairsChannel(params.test_fastqc_reads)
    out = fastqc(reads)
    out.view { i -> "$i" }
}
