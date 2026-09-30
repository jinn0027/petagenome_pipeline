#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_GG2_16S_MAX_MEMORY = 64
def QIIME2_GG2_16S_MAX_THREADS = 32

params.qiime2_gg2_16s_memory = Math.min(params.memory as Integer, QIIME2_GG2_16S_MAX_MEMORY)
params.qiime2_gg2_16s_threads = Math.min(params.threads as Integer, QIIME2_GG2_16S_MAX_THREADS)

// 3. ターゲット領域の指定（デフォルトは V4）
params.qiime2_gg2_16s_target_region = 'v4'

// 4. rrnDBのコピー数探索モード ('right': 右側から細かい階層へ遡る, 'genus': 属レベル固定)
params.qiime2_gg2_rrndb_mode = 'right'

// 5. 機能アノテーション関連パラメータ
if (!params.containsKey('annotation_table')) {
    params.annotation_table = "${params.petagenomeDir}/data/greengenes2/genome_annotations_table.tsv"
}
params.functional_annotations = "KO,MetaCyc"

// リファレンスファイルのパスのデフォルト設定
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_gg2_16s_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_16s_taxonomy = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_gg2_16s_rrndb_stats = "${params.petagenomeDir}/data/rrnDB/rrnDB-5.10_pantaxa_stats_RDP.tsv.zip"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_greengenes2_16s {
    tag "${pair_id} (${target_region})"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_gg2_16s_memory}"
    def threads = "${params.qiime2_gg2_16s_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(table), path(rep_seqs)
        val target_region
        path backbone_fna
        path taxonomy
        path rrndb_stats
        val rrndb_mode
        path annotation_table
        val target_annots

    output:
        tuple val(pair_id), 
              path("${pair_id}/feature-table.tsv"), 
              path("${pair_id}/representatives.fasta"),
              path("${pair_id}/taxonomy.tsv"),
              path("${pair_id}/taxonomy_counts.tsv"),
              path("${pair_id}/*_functional_counts.tsv")

    script:
        def gg2_command = "non-v4-16s" 

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

        echo "Target region: ${target_region}"
        echo "rrnDB search mode: ${rrndb_mode}"

        # 1. Greengenes2 実行（16S用）
        qiime greengenes2 ${gg2_command} \
            --i-table ${table} \
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

        # 5. feature-table.tsv と taxonomy.tsv を結合し、rrnDBコピー数補正 & 全体和正規化を行い、カウント数順（降順）でソートする (タクソノミ)
        python3 ${params.petagenomeDir}/scripts/Python/parse_taxonomy.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${pair_id}/taxonomy_counts.tsv \
            --mode ${rrndb_mode}

        # 6. アノテーションテーブルを紐づけて、rrnDB補正後の存在量を各機能（KO, MetaCycなど）に分配・集計する
        python3 ${params.petagenomeDir}/scripts/Python/parse_functional_profiles.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${annotation_table} \
            "${target_annots}" \
            "${pair_id}" \
            --mode ${rrndb_mode}
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow QIIME2_GREENGENES2_16S_SUB {
    take:
    p
    dada2_out
    target_region
    backbone_fna
    taxonomy
    rrndb_stats
    rrndb_mode
    annotation_table
    target_annots

    main:
    in_ch = dada2_out.map { pair_id, table, rep_seqs, stats, trans ->
        tuple(pair_id, table, rep_seqs)
    }

    out = qiime2_greengenes2_16s(
        p.combine(in_ch).map { p_val, pair_id, table, rep_seqs -> tuple(p_val, pair_id, table, rep_seqs) },
        target_region,
        backbone_fna,
        taxonomy,
        rrndb_stats,
        rrndb_mode,
        annotation_table,
        target_annots
    )

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン用エントリーポイント
// ==========================================
workflow QIIME2_GREENGENES2_16S_ALL {
    p = createNullParamsChannel()
    region_ch = Channel.value(params.qiime2_gg2_16s_target_region)
    mode_ch = Channel.value(params.qiime2_gg2_rrndb_mode)
    
    backbone_ch = Channel.value(file(params.qiime2_gg2_16s_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_16s_taxonomy, checkIfExists: true))
    rrndb_ch = Channel.value(file(params.qiime2_gg2_16s_rrndb_stats, checkIfExists: true))
    annotation_ch = Channel.value(file(params.annotation_table, checkIfExists: true))
    annots_ch = Channel.value(params.functional_annotations)

    dada2_dummy_ch = Channel.fromPath("${params.output}/qiime2_dada2/*/table.qza")
        .map { table_path ->
            def pair_id = table_path.parent.name
            def rep_seqs = file("${table_path.parent}/rep-seqs.qza", checkIfExists: true)
            return tuple(pair_id, table_path, rep_seqs, file('dummy_stats.qza'), file('dummy_trans.qza'))
        }

    QIIME2_GREENGENES2_16S_SUB(
        p,
        dada2_dummy_ch,
        region_ch,
        backbone_ch,
        taxonomy_ch,
        rrndb_ch,
        mode_ch,
        annotation_ch,
        annots_ch
    )
}

workflow {
    QIIME2_GREENGENES2_16S_ALL()
}