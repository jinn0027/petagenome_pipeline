#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.cdhit_cdhit_est_memory = params.memory
params.cdhit_cdhit_est_threads = params.threads

params.cdhit_identity_threshold = 0.95
params.cdhit_global_sequence_identity = 1
params.cdhit_description_length = 150
params.cdhit_word_length = 5
params.cdhit_mask = "NX"

process cdhit_est {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/cdhit/cdhit.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.cdhit_cdhit_est_memory}"
    def threads = "${params.cdhit_cdhit_est_threads}"
    memory "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(read, arity: '1')
    output:
        tuple val(id), path("${id}/out.fasta"), path("${id}/out.fasta.clstr")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        cd-hit-est \\
            -M "${params.cdhit_cdhit_est_memory}000" \\
            -T ${threads} \\
            -c ${params.cdhit_identity_threshold} \\
            -G ${params.cdhit_global_sequence_identity} \\
            -d ${params.cdhit_description_length} \\
            -n ${params.cdhit_word_length} \\
            -mask ${params.cdhit_mask} \\
            -i ${read} \\
            -o ${id}/out.fasta
        """
}

workflow {
    read = createSeqsChannel(params.test_cdhit_read)
    out = cdhit_est(read)
    out.view { i -> "$i" }
}
