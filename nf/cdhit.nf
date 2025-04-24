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
        tuple val(read_id), path(read, arity: '1')

    output:
        tuple val(read_id), path("out/*.fasta", arity: '1')
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
            -o out/${read_id}.fasta
        """
}

workflow {
    read = channel.fromPath(params.test_cdhit_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    out = cdhit_est(read)
    out.view { i -> "$i" }
}
