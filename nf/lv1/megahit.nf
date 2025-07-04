#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process megahit {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/megahit/megahit.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out.contigs.fa", arity: '1')
    script:
        """
        megahit \\
            -m ${params.memory} \\
            -t ${params.threads} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${pair_id} \\
            --out-prefix out
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_megahit_reads, checkIfExists: true)
    out = megahit(reads)
    out.view { i -> "$i" }
}
