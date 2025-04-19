#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.test_spades_separate = true

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}/corrected/*.cor.fastq.gz"), \
	      path("${pair_id}/corrected/unpaired/*.cor.fastq.gz")
    script:
    """
    spades.py \\
        --only-error-correction \\
        --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} \\
        --threads ${params.threads} --memory ${params.memory} -o ${pair_id}
    mkdir -p ${pair_id}/corrected/unpaired
    mv ${pair_id}/corrected/*unpaired*.cor.fastq.gz ${pair_id}/corrected/unpaired
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
   reads = channel.fromFilePairs(params.test_spades_reads, checkIfExists: true)
   if (params.test_spades_separate) {
       error_correction = spades_error_correction(reads)
           .map { pair_id, reads, unpaired -> tuple( pair_id, reads ) }
       //error_correction.view { i -> "$i" }
       spades_assembler = spades_assembler(error_correction)
       //spades_assembler.view { i -> "$i" }
   } else {
       spades = spades(reads)
       //spades.view { i -> "$i" }
   }
}
