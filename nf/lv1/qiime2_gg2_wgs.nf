#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_GG2_WGS_MAX_MEMORY = 64
def QIIME2_GG2_WGS_MAX_THREADS = 32

params.qiime2_gg2_wgs_memory = Math.min(params.memory as Integer, QIIME2_GG2_WGS_MAX_MEMORY)
params.qiime2_gg2_wgs_threads = Math.min(params.threads as Integer, QIIME2_GG2_WGS_MAX_THREADS)

// 3. rrnDBのコピー数探索モード ('right': 右側から細かい階層へ遡る, 'genus': 属レベル固定)
params.qiime2_gg2_rrndb_mode = 'right'

// リファレンスファイルのパスのデフォルト設定
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_gg2_wgs_rrndb_stats = "${params.petagenomeDir}/data/rrnDB/rrnDB-5.10_pantaxa_stats_RDP.tsv.gz"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_greengenes2_wgs {
    tag "${pair_id} (WGS)"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_gg2_wgs_memory}"
    def threads = "${params.qiime2_gg2_wgs_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(rep_seqs)
        path backbone_fna
        path taxonomy
        path rrndb_stats
        val rrndb_mode

    output:
        tuple val(pair_id), 
              path("${pair_id}/feature-table.tsv"), 
              path("${pair_id}/representatives.fasta"),
              path("${pair_id}/taxonomy.tsv"),
              path("${pair_id}/taxonomy_counts.tsv")

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

        echo "rrnDB search mode: ${rrndb_mode}"

        # 1. WGS（Shotgun）向け Greengenes2 実行
        qiime greengenes2 shotgun \
            --i-sequences ${rep_seqs} \
            --i-backbone ${backbone_fna} \
            --o-mapped-table mapped_table.qza \
            --o-representatives representatives.qza \
            --p-threads ${threads}

        # 2. フィーチャーテーブルをTSVに変換
        qiime tools export \
            --input-path mapped_table.qza \
            --output-path exported_table
        
        biom convert \
            -i exported_table/feature-table.biom \
            -o ${pair_id}/feature-table.tsv \
            --to-tsv

        # 3. 代表配列をFASTAに変換
        qiime tools export \
            --input-path representatives.qza \
            --output-path ${pair_id}
        
        mv ${pair_id}/dna-sequences.fasta ${pair_id}/representatives.fasta

        # 4. タクソノミ情報をTSVに変換
        qiime tools export \
            --input-path ${taxonomy} \
            --output-path exported_taxonomy

        if [ -f exported_taxonomy/taxonomy.tsv ]; then
            cp exported_taxonomy/taxonomy.tsv ${pair_id}/taxonomy.tsv
        elif [ -f exported_taxonomy/consensus_assignments.tsv ]; then
            cp exported_taxonomy/consensus_assignments.tsv ${pair_id}/taxonomy.tsv
        else
            find exported_taxonomy -name "*.tsv" -exec cp {} ${pair_id}/taxonomy.tsv
        fi

        # 5. feature-table.tsv と taxonomy.tsv を結合し、rrnDBコピー数補正 & 全体和正規化を行い、カウント数順（降順）でソートする
        python3 ${params.petagenomeDir}/scripts/Python/parse_taxonomy.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${pair_id}/taxonomy_counts.tsv \
            --mode ${rrndb_mode}
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    input_ch // tuple val(ref_id), val(pair_id), path(rep_seqs)
    backbone_fna
    taxonomy
    rrndb_stats
    rrndb_mode

    main:
    // 修正: p.combine(input_ch) で要素数が 4個 (p_val, ref_id, pair_id, rep_seqs) になるため、引数4つで受けてタスク用に再構築
    out = qiime2_greengenes2_wgs(
        p.combine(input_ch).map { p_val, ref_id, pair_id, rep_seqs -> tuple(p_val, pair_id, rep_seqs) },
        backbone_fna,
        taxonomy,
        rrndb_stats,
        rrndb_mode
    )

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン用エントリーポイント
// ==========================================
workflow QIIME2_GREENGENES2_WGS_ALL {
    p = createNullParamsChannel()
    mode_ch = Channel.value(params.qiime2_gg2_rrndb_mode)
    
    backbone_ch = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))
    rrndb_ch = Channel.value(file(params.qiime2_gg2_wgs_rrndb_stats, checkIfExists: true))

    wgs_input_ch = Channel.fromPath("${params.output}/wgs_rep_seqs/*/rep-seqs.qza")
        .map { rep_path ->
            def ref_id = "gg2_backbone"
            def pair_id = rep_path.parent.name
            return tuple(ref_id, pair_id, rep_path)
        }

    QIIME2_GREENGENES2_WGS_SUB(
        p,
        wgs_input_ch,
        backbone_ch,
        taxonomy_ch,
        rrndb_ch,
        mode_ch
    )
}

workflow {
    QIIME2_GREENGENES2_WGS_ALL()
}