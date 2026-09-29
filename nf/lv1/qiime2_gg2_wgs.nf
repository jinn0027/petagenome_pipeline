#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_GG2_WGS_MAX_MEMORY  = 64
def QIIME2_GG2_WGS_MAX_THREADS = 32

params.qiime2_gg2_wgs_memory  = Math.min(params.memory as Integer, QIIME2_GG2_WGS_MAX_MEMORY)
params.qiime2_gg2_wgs_threads = Math.min(params.threads as Integer, QIIME2_GG2_WGS_MAX_THREADS)

// 必須パラメータのチェック
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy     = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_greengenes2_wgs {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_gg2_wgs_memory}"
    def threads = "${params.qiime2_gg2_wgs_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(table)
        path backbone_fna
        path taxonomy

    output:
        tuple val(pair_id), 
              path("${pair_id}/feature-table.tsv"), 
              path("${pair_id}/taxonomy.tsv")

    script:
        """
        # Pythonの非推奨警告を抑制
        export PYTHONWARNINGS="ignore"

        # コンテナ特有のキャッシュ・ホームディレクトリ競合を防ぐ環境変数
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # 1. Greengenes2 を用いたWGSフィーチャーテーブルのフィルタリング・マッピング
        qiime greengenes2 filter-features \
            --i-table ${table} \
            --i-backbone ${backbone_fna} \
            --o-filtered-table filtered_table.qza

        # 2. フィルタリングされたテーブルからGreengenes2ベースのタクソノミーを抽出
        qiime greengenes2 taxonomy-from-table \
            --i-table filtered_table.qza \
            --i-reference-taxonomy ${taxonomy} \
            --o-classification taxonomy.qza

        # 3. フィーチャーテーブルをTSVに変換
        qiime tools export \
            --input-path filtered_table.qza \
            --output-path exported_table
        
        biom convert \
            -i exported_table/feature-table.biom \
            -o ${pair_id}/feature-table.tsv \
            --to-tsv

        # 4. タクソノミー情報をTSVに変換
        qiime tools export \
            --input-path taxonomy.qza \
            --output-path exported_taxonomy

        if [ -f exported_taxonomy/taxonomy.tsv ]; then
            cp exported_taxonomy/taxonomy.tsv ${pair_id}/taxonomy.tsv
        elif [ -f exported_taxonomy/consensus_assignments.tsv ]; then
            cp exported_taxonomy/consensus_assignments.tsv ${pair_id}/taxonomy.tsv
        else
            find exported_taxonomy -name "*.tsv" -exec cp {} ${pair_id}/taxonomy.tsv
        fi
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    woltka_out 
    backbone_fna
    taxonomy

    main:
    in_ch = woltka_out.map { pair_id, table ->
        tuple(pair_id, table)
    }

    out = qiime2_greengenes2_wgs(
        p.combine(in_ch).map { p_val, pair_id, table -> tuple(p_val, pair_id, table) },
        backbone_fna,
        taxonomy
    )

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン用エントリーポイント
// ==========================================
workflow QIIME2_GREENGENES2_WGS_ALL {
    p           = createNullParamsChannel()

    // Nextflowの組み込み機能（checkIfExists: true）で安全にファイルを指定
    backbone_ch = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))

    woltka_dummy_ch = Channel.fromPath("${params.output}/woltka_table/*/table.qza")
        .map { table_path ->
            def pair_id = table_path.parent.name
            return tuple(pair_id, table_path)
        }

    QIIME2_GREENGENES2_WGS_SUB(
        p,
        woltka_dummy_ch,
        backbone_ch,
        taxonomy_ch
    )
}