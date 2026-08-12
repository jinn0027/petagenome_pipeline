#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// 2. MEGAHIT 固有の上限値定義
def MEGAHIT_MAX_MEMORY  = 128 // GB (大規模メタゲノム Graph 構築用の高メモリ上限)
def MEGAHIT_MAX_THREADS = 16  // アセンブリスケール効率とメモリバス競合のバランス上限

// 3. 上限値による動的クリッピング
params.megahit_megahit_memory  = Math.min((params.megahit_megahit_memory  ?: params.memory) as Integer, MEGAHIT_MAX_MEMORY)
params.megahit_megahit_threads = Math.min((params.megahit_megahit_threads ?: params.threads) as Integer, MEGAHIT_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process megahit {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/megahit/megahit.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
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
