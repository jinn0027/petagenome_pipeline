#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cdhit_identity_threshold = 0.95
params.cdhit_global_sequence_identity = 1
params.cdhit_description_length = 150
params.cdhit_word_length = 5
params.cdhit_mask = "NX"

process cdhit_est {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/cdhit/cdhit.sif"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(read, arity: '1')
    output:
        tuple val(id), path("out.fasta"), path("out.fasta.clstr")
    script:
        """
        mkdir -p out
        cd-hit-est \\
            -M "${params.memory}000" \\
            -T ${params.threads} \\
            -c ${params.cdhit_identity_threshold} \\
            -G ${params.cdhit_global_sequence_identity} \\
            -d ${params.cdhit_description_length} \\
            -n ${params.cdhit_word_length} \\
            -mask ${params.cdhit_mask} \\
            -i ${read} \\
            -o out.fasta
        """
}

workflow {
    read = channel.fromPath(params.test_cdhit_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    out = cdhit_est(read)
    out.view { i -> "$i" }
}
