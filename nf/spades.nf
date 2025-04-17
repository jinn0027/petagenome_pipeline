#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.spades_memory = 20

process spades {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/03_spades/spades", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/contigs.fasta")
    script:
    """
    spades.py --meta \\
    --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
    --threads ${params.threads} --memory ${params.spades_memory} -o ${pair_id}
    """
}

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/03_spades/error_correction", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/**.gz")
    script:
    """
    spades.py --only-error-correction \\
    --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
    --threads ${params.threads} --memory ${params.spades_memory} -o ${pair_id}
    """
}

process spades_assembler {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/03_spades/assembly", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/contigs.fasta")
    script:
    """
    spades.py --only-assembler --meta \\
    --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
    --threads ${params.threads} --memory ${params.spades_memory} -o ${pair_id}
    """
}

