#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"
include { error_correction } from "${params.petagenomeDir}/nf/lv2/error_correction"
include { assembly } from "${params.petagenomeDir}/nf/lv2/assembly"
include { pool_contigs } from "${params.petagenomeDir}/nf/lv2/pool_contigs"

workflow bacteriome_pipeline {
  take:
    reads

  main:
    def l_thre = 5000

    fp = fastp(reads)
    ec = error_correction(fp)
    as = assembly(ec.ec.map { pair_id, paired, unpaired -> tuple( pair_id, paired ) }, l_thre)
    pc = pool_contigs(as.flt, l_thre)

  emit:
    fp
    ec.ec
    ec.fqc
    ec.len
    as.flt
    as.len
    as.sts
    as.blstdb
    pc.merged
    pc.clust
    pc.flt
    pc.name
    pc.len
    pc.sts
    pc.blstdb
}

workflow {
  reads = channel.fromFilePairs(params.test_bacteriome_pipeline_reads, checkIfExists: true)
  out = bacteriome_pipeline(reads)
}
