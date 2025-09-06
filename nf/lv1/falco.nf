#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.falco_falco_memory = params.memory
params.falco_falco_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process falco {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/falco/falco.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.falco_falco_memory}"
    def threads = "${params.falco_falco_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        val(p)
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${pair_id}
        falco \\
            -t ${threads} \\
            -o ${pair_id} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.test_falco_reads)
    out = falco(p, reads)
    out.view { i -> "$i" }
}
