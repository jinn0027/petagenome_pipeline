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
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # seqkit pairを実行
        # -1 / -2 で入力を指定
        # -O で出力ディレクトリを指定
        # -f （必要に応じて既存ディレクトリの上書き）
        # -j はグローバルオプションとして渡すか、そのまま指定可能
        
        seqkit pair \
            -j ${threads} \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            -O ${pair_id}_seqkit_tmp \
            -f \
            --quiet

        # seqkit pair の出力ファイル名構造に合わせてリネーム
        # 通常、入力が pair_1.fastq.gz / pair_2.fastq.gz の場合、
        # 出力側にはペアが一致したものが _1.fastq.gz / _2.fastq.gz として格納される
        
        # 出力ファイル名を安全に捕捉して移動
        mv ${pair_id}_seqkit_tmp/*_1*.fastq.gz ${pair_id}/sync_1.fastq.gz
        mv ${pair_id}_seqkit_tmp/*_2*.fastq.gz ${pair_id}/sync_2.fastq.gz
        
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