#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createNullParamsChannel; getParam } from "${params.petagenomeDir}/nf/common/utils"
include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"
include { error_correction } from "${params.petagenomeDir}/nf/lv2/error_correction"
include { assembly } from "${params.petagenomeDir}/nf/lv2/assembly"
include { pool_contigs } from "${params.petagenomeDir}/nf/lv2/pool_contigs"
include { circular_contigs } from "${params.petagenomeDir}/nf/lv2/circular_contigs"
//include { circular_contigs } from "${params.petagenomeDir}/nf/lv2/circular_contigs2"


workflow bacteriome_pipeline {
  take:
    p
    reads
    l_thre
  main:
    fp = fastp(p.combine(reads))
    ec = error_correction(p, fp)
    as = assembly(p, ec.ec.map { pair_id, paired, unpaired -> tuple( pair_id, paired ) }, l_thre)

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

    pc = pool_contigs(p, flt_all, l_thre)
    cc = circular_contigs(p, pc.flt.map { id, fa, name -> [id, fa] })
  emit:
    reads
    fp
    ec.ec
    ec.fqc
    ec.len
    as.flt
    /*
    as.len
    as.sts
    as.blstdb
    */
    pc.merged
    pc.clust
    pc.flt
    pc.name
    pc.len
    pc.sts
    pc.blstdb
    cc.dedupl
}

workflow {
    p = createNullParamsChannel()
    def reads_list = params.test_bacteriome_pipeline_reads.split(';')

    def individual_channels = []
    reads_list.each { reads ->
        def ch = channel.fromFilePairs(reads, checkIfExists: true)
        individual_channels << ch
    }

    def reads_mixed = individual_channels.first()
    individual_channels.tail().each {
        ch -> reads_mixed = reads_mixed.mix(ch)
    }

    index = 0
    def reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index)}_${id}"
        index += 1
        return tuple(new_id, pair)
    }

    out = bacteriome_pipeline(p, reads, params.test_bacteriome_pipeline_lthre)
}
