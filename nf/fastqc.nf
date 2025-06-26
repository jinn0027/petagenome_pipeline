#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqc {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastqc/fastqc.sif"
    publishDir "${params.output}/${task.process}/${pair_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("out")
    script:
        """
        mkdir -p out
        fastqc \\
            -o out\\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_fastqc_reads, checkIfExists: true)
    out = fastqc(reads)
    out.view { i -> "$i" }
}
