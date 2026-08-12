#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// 2. CD-HIT-EST 固有の上限値定義
def CDHIT_EST_MAX_MEMORY  = 64 // GB (大量の代表配列インデックス保持用)
def CDHIT_EST_MAX_THREADS = 16 // OpenMP並列化の競合・飽和点を考慮

// 3. 上限値による動的クリッピング
params.cdhit_cdhit_est_memory  = Math.min((params.cdhit_cdhit_est_memory  ?: params.memory) as Integer, CDHIT_EST_MAX_MEMORY)
params.cdhit_cdhit_est_threads = Math.min((params.cdhit_cdhit_est_threads ?: params.threads) as Integer, CDHIT_EST_MAX_THREADS)

params.cdhit_identity_threshold = 0.95
params.cdhit_global_sequence_identity = 1
params.cdhit_description_length = 150
params.cdhit_word_length = 5
params.cdhit_mask = "NX"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process cdhit_est {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/cdhit/cdhit.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.cdhit_cdhit_est_memory}"
    def threads = "${params.cdhit_cdhit_est_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(id), path(read, arity: '1')
    output:
        tuple val(id), path("${id}/out.fasta"), path("${id}/out.fasta.clstr")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${id}
        cd-hit-est \\
            -M "${gb}000" \\
            -T ${threads} \\
            -c ${getParam(p, 'cdhit_identity_threshold')} \\
            -G ${getParam(p, 'cdhit_global_sequence_identity')} \\
            -d ${getParam(p, 'cdhit_description_length')} \\
            -n ${getParam(p, 'cdhit_word_length')} \\
            -mask ${getParam(p, 'cdhit_mask')} \\
            -i ${read} \\
            -o ${id}/out.fasta
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// CD-HIT-EST 処理の本体
workflow CDHIT_EST_SUB {
    take:
    p
    read

    main:
    // [ p_val, seq_id, seq_path ] にフラット化
    in_ch = p.combine(read).map { p_val, seq_id, seq_path ->
        tuple(p_val, seq_id, seq_path)
    }

    out = cdhit_est(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow CDHIT_EST_ALL {
    p    = createNullParamsChannel()
    read = createSeqsChannel(params.cdhit_read)

    out_ch = CDHIT_EST_SUB(p, read)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    CDHIT_EST_ALL()
}
