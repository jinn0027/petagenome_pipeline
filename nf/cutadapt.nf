#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cutadapt_fwd = "AATGATACGGCGACCACCGAGAUCTACAC"
params.cutadapt_rev = "CAAGCAGAAGACGGCATACGAGAT"
params.cutadapt_minimum_length = 50

process cutadapt {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/cutadapt/cutadapt.sif"
    publishDir "${params.output}/cutadapt/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("${pair_id}_*.fastq.gz", arity: '2')
    script:
        """
        cutadapt \\
            -a ${params.cutadapt_fwd} \\
            -g ${params.cutadapt_rev} \\
            -o ${pair_id}_cutadapt_out1.fastq.gz \\
            -p ${pair_id}_cutadapt_out2.fastq.gz \\
            --minimum-length ${params.cutadapt_minimum_length} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
   reads = channel.fromFilePairs(params.test_cutadapt_reads, checkIfExists: true)
   out = cutadapt(reads)
   //out.view { i -> "$i" }
}
