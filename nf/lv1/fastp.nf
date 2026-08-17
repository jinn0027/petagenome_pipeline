#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. このモジュール・タスク固有の推奨・上限値（ローカル定数として定義）
def FASTP_MAX_MEMORY  = 16  // fastp はメモリをほぼ使わないため 16GB 上限
def FASTP_MAX_THREADS = 8   // fastp のスレッド数は 8 以上で並列効率が飽和する

// 3. 上限値による動的クリッピング
params.fastp_fastp_memory  = Math.min(params.memory as Integer, FASTP_MAX_MEMORY)
params.fastp_fastp_threads = Math.min(params.threads as Integer, FASTP_MAX_THREADS)

params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process fastp {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/fastp/fastp.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.fastp_fastp_memory}"
    def threads = "${params.fastp_fastp_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}/out_{1,2}.fastq*", arity: '2')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        fastp \
            -w ${threads} \
            --detect_adapter_for_pe \
            --low_complexity_filter \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${pair_id}/out_1.fastq.gz \
            -O ${pair_id}/out_2.fastq.gz \
            --cut_front --cut_tail \
            --cut_mean_quality ${getParam(p, params, 'fastp_cut_mean_quality')} \
            --length_required ${getParam(p, params, 'fastp_reads_minlength')} \
            --json ${pair_id}_fastp.json \
            --html ${pair_id}_fastp.html
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// fastp 処理の本体
workflow FASTP_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = fastp(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow FASTP_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.fastp_reads)

    out_ch = FASTP_SUB(p, reads)
    out_ch.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    FASTP_ALL()
}