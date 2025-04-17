#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${projectDir}/fastp.nf"
include { spades } from "${projectDir}/spades.nf"

params.reads = "${projectDir}/../modules/test/s_6_{1,2}.fastq.gz"
params.output = "output"
//params.threads = 4
//params.fastp_cut_mean_quality = 15
//params.fastp_reads_minlength = 15
//params.fastp_memory = 10
////params.spades_memory = 20

workflow {
   raw_short_reads = channel.fromFilePairs(params.reads, checkIfExists: true)
   //raw_short_reads.view{ raw_short_pairs -> "$pairs"}
   fastp = fastp(raw_short_reads)
   spades = spades(fastp)
}
