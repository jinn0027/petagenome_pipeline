#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"
include { spades_e2e; spades_error_correction; spades_assembler } from "${params.petagenomeDir}/nf/lv1/spades"

workflow {
   raw_short_reads = channel.fromFilePairs(params.test_main_reads, checkIfExists: true)
   fastp = fastp(raw_short_reads)
   //spades = spades_error_correction(fastp)
   spades = spades_assembler(fastp)
}
