#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createNullParamsChannel; getParam } from "${params.petagenomeDir}/nf/common/utils"
include { fastp } from "${params.petagenomeDir}/nf/lv1/fastp"

workflow bacteriome_pipeline {
  take:
    p
    reads
  main:
    fp = fastp(p.combine(reads))
  emit:
    reads
    fp
}

workflow {
    p = createNullParamsChannel()
    def reads_list = params.test_bacteriome_pipeline_reads.split(';')

    // 1. パスごとにチャンネルを作ってリストに格納
    def individual_channels = []
    reads_list.each { reads ->
        def ch = channel.fromFilePairs(reads, checkIfExists: true)
        individual_channels << ch
    }

    // 2. 複数のチャンネルを1つに統合
    def reads_mixed = individual_channels.first()
    individual_channels.tail().each { ch ->
        reads_mixed = reads_mixed.mix(ch)
    }

    // 3. IDに連番を付与
    def index = 0
    def reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index++)}_${id}"
        return tuple(new_id, pair)
    }

    out = bacteriome_pipeline(p, reads)
}
