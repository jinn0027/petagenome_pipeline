#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. このモジュール・タスク固有の推奨・上限値
def SYNC_MAX_MEMORY  = 32  
def SYNC_MAX_THREADS = 16  // seqkit はマルチスレッド処理に対応している

// 3. 上限値による動的クリッピング
params.sync_pair_memory  = Math.min(params.memory as Integer, SYNC_MAX_MEMORY)
params.sync_pair_threads = Math.min(params.threads as Integer, SYNC_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process sync_pair {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    
    def gb      = "${params.sync_pair_memory}"
    def threads = "${params.sync_pair_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus   params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("${pair_id}/sync_{1,2}.fastq.gz", arity: '2')

    script:
        // 入力パスからファイル名（例: ERR1620255_..._R1_001.fastq.gz）を抽出しておく
        def r1_name = new File("${reads[0]}").name
        def r2_name = new File("${reads[1]}").name
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        seqkit pair \
            -j ${threads} \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            -O ${pair_id}_seqkit_tmp \
            -f \
            --quiet

        # seqkit pair が出力した元のファイル名でキャプチャしてリネーム移動
        mv ${pair_id}_seqkit_tmp/${r1_name} ${pair_id}/sync_1.fastq.gz
        mv ${pair_id}_seqkit_tmp/${r2_name} ${pair_id}/sync_2.fastq.gz
        
        # 一時ディレクトリの削除
        rm -rf ${pair_id}_seqkit_tmp
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow SYNC_PAIR_SUB {
    take:
    p
    reads

    main:
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = sync_pair(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. エントリーポイント (-entry)
// ==========================================
workflow SYNC_PAIR_ALL {
    p      = createNullParamsChannel()
    reads  = createPairsChannel(params.sync_pair_reads)

    out_ch = SYNC_PAIR_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

workflow {
    SYNC_PAIR_ALL()
}