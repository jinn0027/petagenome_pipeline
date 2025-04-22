#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cdhit_identity_threshold = 0.95
params.cdhit_global_sequence_identity = 1
params.cdhit_description_length = 150
params.cdhit_word_length = 5
params.cdhit_mask = "NX"

process cdhit_est {
    tag "${read_id}"
    container = "${params.petagenomeDir}/modules/cdhit/cdhit.sif"
    publishDir "${params.output}/cdhit/${read_id}", mode: 'copy'
    input:
        tuple val(read_id), path(read)

    output:
        tuple val(read_id), path("${read_id}_cdhit_out.fasta")
    script:
    """
    cd-hit-est \\
        -M ${params.memory} \\
        -T ${params.threads} \\
        -c ${params.cdhit_identity_threshold} \\
        -G ${params.cdhit_global_sequence_identity} \\
        -d ${params.cdhit_description_length} \\
        -n ${params.cdhit_word_length} \\
        -mask ${params.cdhit_mask} \\
        -i ${read} \\
        -o ${read_id}_cdhit_out.fasta
    """
}

workflow {
    read = channel.fromPath(params.test_cdhit_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    cdhit_est = cdhit_est(read)
    //cdhit_est.view { i -> "$i" }
}
