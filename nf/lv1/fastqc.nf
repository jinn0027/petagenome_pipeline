#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { printProcessProfile; createPairsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.fastqc_fastqc_memory = params.memory
params.fastqc_fastqc_threads = params.threads

process fastqc {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastqc/fastqc.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.fastqc_fastqc_memory} GB"
    cpus "${params.fastqc_fastqc_threads}"

    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        printProcessProfile(task)
        """
        export XDG_CACHE_HOME=\$(pwd)/.cache
        mkdir -p .cache
        mkdir -p ${pair_id}
        fastqc \\
            -t ${params.fastqc_fastqc_threads} \\
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
