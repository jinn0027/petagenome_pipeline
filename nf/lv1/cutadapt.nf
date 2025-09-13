#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cutadapt_cutadapt_memory = params.memory
params.cutadapt_cutadapt_threads = params.threads

params.cutadapt_fwd = "AATGATACGGCGACCACCGAGAUCTACAC"
params.cutadapt_rev = "CAAGCAGAAGACGGCATACGAGAT"
params.cutadapt_minimum_length = 50

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process cutadapt {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/cutadapt/cutadapt.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.cutadapt_cutadapt_memory}"
    def threads = "${params.cutadapt_cutadapt_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out_{1,2}.fastq", arity: '2')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        cutadapt \\
            -a ${getParam(p, 'cutadapt_fwd')} \\
            -g ${getParam(p, 'cutadapt_rev')} \\
            -o ${pair_id}/out_1.fastq \\
            -p ${pair_id}/out_2.fastq \\
            --minimum-length ${getParam(p, 'cutadapt_minimum_length')} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

workflow {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.test_cutadapt_reads)
    out = cutadapt(p.combine(reads))
    //out.view { i -> "$i" }
}
