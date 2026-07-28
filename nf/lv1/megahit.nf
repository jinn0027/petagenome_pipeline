#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.megahit_megahit_memory = params.memory
params.megahit_megahit_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process megahit {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/megahit/megahit.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.megahit_megahit_memory}"
    def threads = "${params.megahit_megahit_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out.contigs.fa", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        megahit \\
            -m "${gb}000000000" \\
            -t ${threads} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${pair_id} \\
            --out-prefix out
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// MEGAHIT 処理の本体
workflow MEGAHIT_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = megahit(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow MEGAHIT_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.megahit_reads)

    out_ch = MEGAHIT_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    MEGAHIT_ALL()
}
