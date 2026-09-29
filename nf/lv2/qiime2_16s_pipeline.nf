#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory = 32
params.threads = 8

// パラメータの初期化
params.qiime2_dada2_trim_left_f = 0
params.qiime2_dada2_trim_left_r = 0
params.qiime2_dada2_trunc_len_f = 230
params.qiime2_trunc_len_r = 230

params.qiime2_gg2_16s_target_region = 'v4'
params.qiime2_gg2_16s_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_16s_taxonomy = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"

// 必須パラメータのチェック
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it."
}

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// 下位モジュール（lv1）のインポート（qiime2_gg2_16sに対応）
include { QIIME2_DADA2_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_dada2"
include { QIIME2_GREENGENES2_16S_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_gg2_16s"

// ==========================================
// サブワークフロー（処理の本体）
// ==========================================
workflow QIIME2_16S_PIPELINE_SUB {
    take:
    p
    reads
    target_region
    backbone_fna
    taxonomy

    main:
    // A. DADA2 によるデノイジング
    dada2_out = QIIME2_DADA2_SUB(p, reads)

    // B. Greengenes2 (16S用) による系統配置・タクソノミ付与
    gg2_out = QIIME2_GREENGENES2_16S_SUB(
        p,
        dada2_out,
        target_region,
        backbone_fna,
        taxonomy
    )

    emit:
    gg2_out = gg2_out
}

// ==========================================
// コマンドライン用エントリーポイント (-entry)
// ==========================================
workflow QIIME2_16S_PIPELINE_ALL {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.qiime2_reads)
    region_ch = Channel.value(params.qiime2_gg2_16s_target_region)

    // Nextflowの標準機能（checkIfExists: true）で安全にファイル存在チェックを行う
    backbone_ch = Channel.value(file(params.qiime2_gg2_16s_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_16s_taxonomy, checkIfExists: true))

    out_ch = QIIME2_16S_PIPELINE_SUB(
        p,
        reads,
        region_ch,
        backbone_ch,
        taxonomy_ch
    )

    out_ch.gg2_out.view { i -> "QIIME2 16S PIPELINE OUT: $i" }
}

workflow {
    QIIME2_16S_PIPELINE_ALL()
}