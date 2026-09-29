#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 32
params.threads = 8

// パラメータの初期化
params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy     = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
// Woltka用のマップファイルなどのデフォルトパスも必要に応じて設定
params.qiime2_woltka_wol_map       = "${params.petagenomeDir}/data/greengenes2/wol_map.txt"

// 必須パラメータのチェック
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it."
}

include { createNullParamsChannel; createPairsChannel; createSeqsChannel; apptainerContainerOptions } from "${params.petagenomeDir}/nf/common/utils"

// 下位モジュール（lv1）のインポート
include { FASTP_SUB } from "${params.petagenomeDir}/nf/lv1/fastp"
include { BUILD_REF_DB_SUB; MAP_SUB } from "${params.petagenomeDir}/nf/lv1/bowtie2"
include { QIIME2_WOLTKA_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_woltka"
include { QIIME2_GREENGENES2_WGS_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_gg2_wgs"

// ==========================================
// .qza から FASTA へのエクスポートプロセス
// ==========================================
process export_qiime2_fasta {
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    input:
        path(qza_file)
    output:
        path("backbone_exported.fasta")
    script:
        """
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        qiime tools export \
            --input-path ${qza_file} \
            --output-path exported_ref

        if [ -f exported_ref/dna-sequences.fasta ]; then
            cp exported_ref/dna-sequences.fasta backbone_exported.fasta
        elif ls exported_ref/*.fasta 1> /dev/null 2>&1; then
            cp exported_ref/*.fasta backbone_exported.fasta
        else
            cp exported_ref/*.fna backbone_exported.fasta
        fi
        """
}

// ==========================================
// サブワークフロー（処理の本体）
// ==========================================
workflow QIIME2_WGS_PIPELINE_SUB {
    take:
    p
    reads
    backbone
    taxonomy
    wol_map

    main:

    // A. fastp による品質管理・フィルタリング
    fastp_out = FASTP_SUB(p, reads)

    // B-0. .qza をプレーンな FASTA にエクスポート
    exported_ref = export_qiime2_fasta(backbone)
    // createSeqsChannel が作る構造にあわせて、[ref_id, [path]] の形に明示的に構築する
    exported_ref_ch = exported_ref.map { fasta_path ->
        def ref_id = "gg2_backbone" // 必要に応じた識別子
        return tuple(ref_id, [fasta_path])
    }
    
    // B-1. Bowtie2 リファレンス DB の作成（エクスポートされた FASTA を使用）
    ref_db = BUILD_REF_DB_SUB(p, exported_ref_ch)

    // B-2. Bowtie2 によるリファレンスへのマッピング
    bowtie2_out = MAP_SUB(p, ref_db, fastp_out)

    // C. Woltka によるプロファイリング・カウントテーブル集計
    woltka_out = QIIME2_WOLTKA_SUB(p, bowtie2_out, wol_map, backbone, taxonomy)

    // D. Greengenes2 (WGS用) によるフィルタリング・タクソノミ付与
    gg2_out = QIIME2_GREENGENES2_WGS_SUB(p, woltka_out, backbone, taxonomy)

    emit:
    gg2_out = gg2_out
}

// ==========================================
// コマンドライン用エントリーポイント (-entry)
// ==========================================
workflow QIIME2_WGS_PIPELINE_ALL {
    p           = createNullParamsChannel()
    reads       = createPairsChannel(params.qiime2_reads)
    
    // Nextflowの標準機能（checkIfExists: true）で安全にファイル存在チェックを行う
    backbone_ch = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))
    wol_map_ch  = Channel.value(file(params.qiime2_woltka_wol_map, checkIfExists: true))

    out_ch = QIIME2_WGS_PIPELINE_SUB(
        p,
        reads,
        backbone_ch,
        taxonomy_ch,
        wol_map_ch
    )

    out_ch.gg2_out.view { i -> "QIIME2 WGS PIPELINE OUT: $i" }
}

workflow {
    QIIME2_WGS_PIPELINE_ALL()
}