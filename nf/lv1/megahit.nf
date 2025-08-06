#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createPairsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.megahit_megahit_memory = params.memory
params.megahit_megahit_threads = params.threads

process megahit {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/megahit/megahit.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.megahit_megahit_memory} GB"
    def threads = "${params.megahit_megahit_threads}"
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out.contigs.fa", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        megahit \\
            -m "${params.megahit_megahit_memory}000000000" \\
            -t ${threads} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${pair_id} \\
            --out-prefix out
        """
}

workflow {
    reads = createPairsChannel(params.test_megahit_reads)
    out = megahit(reads)
    out.view { i -> "$i" }
}
