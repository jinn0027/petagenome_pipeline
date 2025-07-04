#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.test_spades_e2e = false

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/corrected/paired/*.cor.fastq", arity: '0..2'), \
              path("${pair_id}/corrected/unpaired/*.cor.fastq", arity: '0..*')
    script:
        """
        mkdir ${pair_id}
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --only-error-correction \\
            --disable-gzip-output \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        mkdir -p ${pair_id}/corrected/paired ${pair_id}/corrected/unpaired 
        mv ${pair_id}/corrected/*unpaired*.cor.fastq ${pair_id}/corrected/unpaired
        mv ${pair_id}/corrected/*.cor.fastq ${pair_id}/corrected/paired
        """
}

process spades_error_correction_gzip_output {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/corrected/paired/*.cor.fastq.gz", arity: '0..2'), \
              path("${pair_id}/corrected/unpaired/*.cor.fastq.gz", arity: '0..*')
    script:
        """
        mkdir ${pair_id}
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --only-error-correction \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        mkdir -p ${pair_id}/corrected/paired ${pair_id}/corrected/unpaired 
        mv ${pair_id}/corrected/*unpaired*.cor.fastq.gz ${pair_id}/corrected/unpaired
        mv ${pair_id}/corrected/*.cor.fastq.gz ${pair_id}/corrected/paired
        """
}

process spades_assembler {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/scaffolds.fasta", arity: '0..*'), \
              path("${pair_id}/contigs.fasta", arity: '0..*')
    script:
        """
        mkdir -p ${pair_id}
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --only-assembler \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        """
}

process spades_e2e {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/scaffolds.fasta", arity: '0..*'), \
              path("${pair_id}/contigs.fasta", arity: '0..*')
    script:
        """
        mkdir -p ${pair_id}
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        """
}

workflow {
   reads = channel.fromFilePairs(params.test_spades_reads, checkIfExists: true)
   if (params.test_spades_e2e) {
       out = spades_e2e(reads)
       out.view { i -> "$i" }
   } else {
       ec = spades_error_correction(reads)
           .map { pair_id, paired, unpaired -> tuple( pair_id, paired ) }
       ec.view { i -> "$i" }
       out = spades_assembler(ec)
       out.view { i -> "$i" }
   }
}
