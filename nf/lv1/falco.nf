#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process falco {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/falco/falco.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        mkdir -p ${pair_id}
        falco \\
            -t ${params.threads} \\
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
