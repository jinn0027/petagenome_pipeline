#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process megahit {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/megahit/megahit.sif"
    publishDir "${params.output}/megahit/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("out/${pair_id}.contigs.fa")
    script:
    """
    megahit \\
        -m ${params.memory} \\
        -t ${params.threads} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o out \\
        --out-prefix ${pair_id}
    """
}

workflow {
    reads = channel.fromFilePairs(params.test_megahit_reads, checkIfExists: true)
    megahit = megahit(reads)
    //megahit.view { i -> "$i" }
}
