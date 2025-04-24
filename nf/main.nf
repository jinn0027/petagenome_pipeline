#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${projectDir}/fastp"
include { spades_e2e; spades_error_correction; spades_assembler } from "${projectDir}/spades"

workflow {
   raw_short_reads = channel.fromFilePairs(params.reads, checkIfExists: true)
   fastp = fastp(raw_short_reads)
   //spades = spades_error_correction(fastp)
   spades = spades_assembler(fastp)
}
