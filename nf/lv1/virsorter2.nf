#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// 2. VirSorter2 固有の上限値定義
def VIRSORTER2_MAX_MEMORY  = 64 // GB (HMM DB・機械学習モデル特徴量抽出用)
def VIRSORTER2_MAX_THREADS = 16 // HMMER・マルチプロセス並列化の I/O スケール上限

// 3. 上限値による動的クリッピング
params.virsorter2_virsorter2_memory  = Math.min(params.memory as Integer, VIRSORTER2_MAX_MEMORY)
params.virsorter2_virsorter2_threads = Math.min(params.threads as Integer, VIRSORTER2_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process virsorter2 {
    tag "${read_id}"
    def local_db = "/opt/VirSorter2/db"
    container = "${params.petagenomeDir}/modules/virsorter2/virsorter2.sif"
    containerOptions = "${params.apptainerRunOptions} -B ${params.virsorter2_db}:${local_db}"
    //containerOptions = "${params.apptainerRunOptions} -B ${params.virsorter2_db}:${local_db} --writable-tmpfs"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.virsorter2_virsorter2_memory}"
    def threads = "${params.virsorter2_virsorter2_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("${read_id}/final-viral-boundary.tsv"), path("${read_id}/final-viral-score.tsv")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        virsorter \\
            config \\
            --init-source \\
            --db-dir=${local_db}
        virsorter \\
            run \\
            -j ${threads} \\
            -w ${read_id} \\
            -i ${read}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// VirSorter2 処理の本体
workflow VIRSORTER2_SUB {
    take:
    p
    read

    main:
    // [ p_val, seq_id, seq_path ] にフラット化
    in_ch = p.combine(read).map { p_val, seq_id, seq_path ->
        tuple(p_val, seq_id, seq_path)
    }

    out = virsorter2(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow VIRSORTER2_ALL {
    p    = createNullParamsChannel()
    read = createSeqsChannel(params.virsorter2_read)

    out_ch = VIRSORTER2_SUB(p, read)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    VIRSORTER2_ALL()
}
