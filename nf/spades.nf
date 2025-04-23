#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.test_spades_e2e = false

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("out/corrected/paired/*.cor.fastq.gz", arity: '2'), \
	      path("out/corrected/unpaired/*.cor.fastq.gz")
    script:
        """
        mkdir out
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --only-error-correction \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o out
        mkdir -p out/corrected/paired out/corrected/unpaired 
        mv out/corrected/*unpaired*.cor.fastq.gz out/corrected/unpaired
        mv out/corrected/*.cor.fastq.gz out/corrected/paired
        """
}

process spades_assembler {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("out/scaffolds.fasta", arity: '0..*'), \
              path("out/contigs.fasta", arity: '0..*')
    script:
        """
        mkdir -p out
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --only-assembler \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o out
        """
}

process spades_e2e {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    publishDir "${params.output}/spades/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("out/scaffolds.fasta", arity: '0..*'), \
              path("out/contigs.fasta", arity: '0..*')
    script:
        """
        mkdir -p out
        spades.py \\
            --threads ${params.threads} \\
            --memory ${params.memory} \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o out
        """
}

workflow {
   reads = channel.fromFilePairs(params.test_spades_reads, checkIfExists: true)
   if (params.test_spades_e2e) {
       out = spades_e2e(reads)
       out.view { i -> "$i" }
   } else {
       ec = spades_error_correction(reads)
           .map { pair_id, reads, unpaired -> tuple( pair_id, reads ) }
       ec.view { i -> "$i" }
       out = spades_assembler(ec)
       out.view { i -> "$i" }
   }
}
