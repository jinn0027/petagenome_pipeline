#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. FastQC 固有の上限値定義
def FASTQC_MAX_MEMORY  = 16 // GB (JVMメモリ＋複数ファイル同時処理用)
def FASTQC_MAX_THREADS = 8  // ファイル並列数・ディスク I/O 競合の限界点

// 3. 上限値による動的クリッピング
params.fastqc_fastqc_memory  = Math.min(params.memory as Integer, FASTQC_MAX_MEMORY)
params.fastqc_fastqc_threads = Math.min(params.threads as Integer, FASTQC_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process fastqc {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastqc/fastqc.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.fastqc_fastqc_memory}"
    def threads = "${params.fastqc_fastqc_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        export XDG_CACHE_HOME=\$(pwd)/.cache
        mkdir -p .cache
        mkdir -p ${pair_id}
        fastqc \\
            -t ${threads} \\
            -o ${pair_id} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// FastQC 処理の本体
workflow FASTQC_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = fastqc(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow FASTQC_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.fastqc_reads)

    out_ch = FASTQC_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    FASTQC_ALL()
}
