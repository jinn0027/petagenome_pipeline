#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.spades_spades_error_correction_memory = params.memory
params.spades_spades_error_correction_threads = params.threads

params.spades_spades_error_correction_gzip_output_memory = params.memory
params.spades_spades_error_correction_gzip_output_threads = params.threads

params.spades_spades_assembler_memory = params.memory
params.spades_spades_assembler_threads = params.threads

params.spades_spades_e2e_memory = params.memory
params.spades_spades_e2e_threads = params.threads

params.spades_error_correction_threads = params.threads
params.spades_error_correction_memory = params.memory
params.test_spades_e2e = false

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.spades_spades_error_correction_memory} GB"
    cpus "${params.spades_spades_error_correction_threads}"

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
            --memory ${params.spades_spades_error_correction_memory} \\
            --threads ${params.spades_spades_error_correction_threads} \\
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
    memory "${params.spades_spades_error_correction_gzip_output_memory} GB"
    cpus "${params.spades_spades_error_correction_gzip_output_threads}"

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
            --memory ${params.spades_spades_error_correction_gzip_output_memory} \\
            --threads ${params.spades_spades_error_correction_gzip_output_threads} \\
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
    memory "${params.spades_spades_assembler_memory} GB"
    cpus "${params.spades_spades_assembler_threads}"

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
            --memory ${params.spades_spades_assembler_memory} \\
            --threads ${params.spades_spades_assembler_threads} \\
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
            --memory ${params.spades_spades_e2e_memory} \\
            --threads ${params.spades_spades_e2e_threads} \\
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
       ec = spades_error_correction_gzip_output(reads)
           .map { pair_id, paired, unpaired -> tuple( pair_id, paired ) }
       ec.view { i -> "$i" }
       out = spades_assembler(ec)
       out.view { i -> "$i" }
   }
}
