#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory = 32
params.threads = 8

// 2. このモジュール・タスク固有の推奨・上限値（ローカル定数として定義）
def QIIME2_MAX_MEMORY = 64  // DADA2はメモリを多く消費するため大きめ
def QIIME2_MAX_THREADS = 32 // スレッド数上限

// 3. 上限値による動的クリッピング
params.qiime2_dada2_memory = Math.min(params.memory as Integer, QIIME2_MAX_MEMORY)
params.qiime2_dada2_threads = Math.min(params.threads as Integer, QIIME2_MAX_THREADS)

// QIIME 2 DADA2 パラメータのデフォルト
params.qiime2_trim_left_f = 0
params.qiime2_trim_left_r = 0
params.qiime2_trunc_len_f = 230
params.qiime2_trunc_len_r = 230

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_dada2 {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.qiime2_dada2_memory}"
    def threads = "${params.qiime2_dada2_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
        
    output:
        tuple val(pair_id), 
              path("${pair_id}/table.qza"), 
              path("${pair_id}/rep-seqs.qza"), 
              path("${pair_id}/denoising-stats.qza"), 
              path("${pair_id}/base-transition-stats.qza")

    script:
        """
        # Pythonの非推奨警告を抑制
        export PYTHONWARNINGS="ignore"

        # コンテナ特有のキャッシュ・ホームディレクトリ競合を防ぐための環境変数退避
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # 1. 作業ディレクトリ（\$PWD）を利用して絶対パスを構築し、TSVとして書き込む
        printf "sample-id\\tforward-absolute-filepath\\treverse-absolute-filepath\\n" > manifest.tsv
        printf "%s\\t\$PWD/%s\\t\$PWD/%s\\n" "${pair_id}" "${reads[0]}" "${reads[1]}" >> manifest.tsv

        # 2. FASTQのインポート (Artifact生成)
        qiime tools import \
            --type SampleData[PairedEndSequencesWithQuality] \
            --input-path manifest.tsv \
            --output-path ${pair_id}/demux.qza \
            --input-format PairedEndFastqManifestPhred33V2

        # 3. DADA2によるデノイジング・ASV生成
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs ${pair_id}/demux.qza \
            --p-trim-left-f ${getParam(p, params, 'qiime2_trim_left_f')} \
            --p-trim-left-r ${getParam(p, params, 'qiime2_trim_left_r')} \
            --p-trunc-len-f ${getParam(p, params, 'qiime2_trunc_len_f')} \
            --p-trunc-len-r ${getParam(p, params, 'qiime2_trunc_len_r')} \
            --o-table ${pair_id}/table.qza \
            --o-representative-sequences ${pair_id}/rep-seqs.qza \
            --o-denoising-stats ${pair_id}/denoising-stats.qza \
            --o-base-transition-stats ${pair_id}/base-transition-stats.qza \
            --p-n-threads ${threads}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================
workflow QIIME2_DADA2_SUB {
    take:
    p
    reads

    main:
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = qiime2_dada2(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================
workflow QIIME2_DADA2_ALL {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.qiime2_reads)

    out_ch = QIIME2_DADA2_SUB(p, reads)
    out_ch.view { i -> "$i" }
}

workflow {
    QIIME2_DADA2_ALL()
}