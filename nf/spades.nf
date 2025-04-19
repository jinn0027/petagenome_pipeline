#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/**.gz")
    script:
    """
    spades.py \\
        --only-error-correction \\
        --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
        --threads ${params.threads} --memory ${params.memory} -o ${pair_id}
    """
}

process spades_assembler {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/scaffolds.fasta"), path("${pair_id}/contigs.fasta")
    script:
    """
    spades.py \\
        --only-assembler --meta \\
        --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
        --threads ${params.threads} --memory ${params.memory} -o ${pair_id}
    """
}

process spades {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/scaffolds.fasta"), path("${pair_id}/contigs.fasta")
    script:
    """
    spades.py \\
        --meta \\
        --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
        --threads ${params.threads} --memory ${params.memory} -o ${pair_id}
    """
}

workflow {
   test_spades_reads = channel.fromFilePairs(params.test_spades_reads, checkIfExists: true)
   spades_error_correction(test_spades_reads)
   spades_assembler(test_spades_reads)
}
