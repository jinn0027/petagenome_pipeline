#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.falco_falco_memory = params.memory
params.falco_falco_threads = params.threads

process falco {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/falco/falco.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.falco_falco_memory} GB"
    cpus "${params.falco_falco_threads}"

    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        mkdir -p ${pair_id}
        falco \\
            -t ${params.falco_falco_threads} \\
            -o ${pair_id} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_falco_reads, checkIfExists: true)
    out = falco(reads)
    out.view { i -> "$i" }
}
