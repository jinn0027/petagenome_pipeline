#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"
include { error_correction } from "${params.petagenomeDir}/nf/lv2/error_correction"
include { assembly } from "${params.petagenomeDir}/nf/lv2/assembly"
include { pool_contigs } from "${params.petagenomeDir}/nf/lv2/pool_contigs"

workflow bacteriome_pipeline {
  take:
    reads
    l_thre

  main:
    fp = fastp(reads)
    ec = error_correction(fp)
    as = assembly(ec.ec.map { pair_id, paired, unpaired -> tuple( pair_id, paired ) }, l_thre)

    def flt_all

    if (true) {
        flt_collected = as.flt.collect(flat: false, sort: true)
        flt_all = flt_collected.map { list_of_tuples ->
            def first_key = list_of_tuples[0][0] // リストの最初のタプルの最初の要素を結果のキーにする
            def all_contigs = []
            list_of_tuples.each { key, contigs ->
                all_contigs << contigs
            }
            [first_key, all_contigs]
        }
    } else {
        flt_all = as.flt.groupTuple()
    }

    pc = pool_contigs(flt_all, l_thre)

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
  def reads_list = params.test_bacteriome_pipeline_reads.split(':')

  def individual_channels = []
  reads_list.each { reads ->
    def ch = channel.fromFilePairs(reads, checkIfExists: true)
    individual_channels << ch
  }

  def reads = individual_channels.first()
  individual_channels.tail().each {
    ch -> reads = reads.mix(ch)
  }

  out = bacteriome_pipeline(reads, params.test_bacteriome_pipeline_lthre)
}
