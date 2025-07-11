#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"
include { error_correction } from "${params.petagenomeDir}/nf/lv2/error_correction"
include { assembly } from "${params.petagenomeDir}/nf/lv2/assembly"

workflow {
    reads = channel.fromFilePairs(params.test_bacteriome_pipeline_reads, checkIfExists: true)
    fastp = fastp(reads)
    //fastp.view { i -> "$i" }
    err_corr = error_correction(fastp)
    //err_corr.ec.view{ i -> "$i" }
    //err_corr.fqc.view{ i -> "$i" }
    //err_corr.len.view{ i -> "$i" }
    asm = assembly(err_corr.ec)
    asm.asm.view{ i -> "$i" }
    asm.flt.view{ i -> "$i" }
    asm.len.view{ i -> "$i" }
    asm.sts.view{ i -> "$i" }
    asm.blstdb.view{ i -> "$i" }
}
