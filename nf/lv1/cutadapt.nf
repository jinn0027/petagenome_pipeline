#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createPairsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.cutadapt_cutadapt_memory = params.memory
params.cutadapt_cutadapt_threads = params.threads

params.cutadapt_fwd = "AATGATACGGCGACCACCGAGAUCTACAC"
params.cutadapt_rev = "CAAGCAGAAGACGGCATACGAGAT"
params.cutadapt_minimum_length = 50

process cutadapt {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/cutadapt/cutadapt.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.cutadapt_cutadapt_memory} GB"
    cpus "${params.cutadapt_cutadapt_threads}"

    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out_{1,2}.fastq", arity: '2')
    script:
        """
        mkdir -p ${pair_id}
        cutadapt \\
            -a ${params.cutadapt_fwd} \\
            -g ${params.cutadapt_rev} \\
            -o ${pair_id}/out_1.fastq \\
            -p ${pair_id}/out_2.fastq \\
            --minimum-length ${params.cutadapt_minimum_length} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    reads = createPairsChannel(params.test_cutadapt_reads)
    out = cutadapt(reads)
    //out.view { i -> "$i" }
}
