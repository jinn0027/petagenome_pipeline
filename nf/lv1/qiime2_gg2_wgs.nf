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

// 4. 機能アノテーション関連パラメータ
if (!params.containsKey('annotation_table')) {
    params.annotation_table = "${params.petagenomeDir}/data/greengenes2/genome_annotations_table.tsv"
}
params.functional_annotations = "KO,MetaCyc"

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
        tuple val(p), val(pair_id), path(sam_file)
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

        # 1. Bowtie2のSAMファイルからリファレンス配列ごとのカウントを集計し、feature-table.tsv (TSVフォーマットのBIOMライク形式) を作成する
        python3 -c "
        import sys
        from collections import Counter

        counts = Counter()
        sam_path = '${sam_file}'

        with open(sam_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) > 2:
                    ref_id = parts[2]
                    if ref_id != '*':
                        counts[ref_id] += 1

        # feature-table.tsv を出力 (QIIME形式のタブ区切りテーブル：行名=FeatureID, 列=サンプル名)
        sample_id = '${pair_id}'
        with open('${pair_id}/feature-table.tsv', 'w') as out:
            out.write(f'# Constructed from biom file\\n#OTU ID\\t{sample_id}\\n')
            for ref_id, count in counts.items():
                out.write(f'{ref_id}\\t{count}\\n')
        "

        # 2. 代表配列 (.qza) から代表配列 FASTA をエクスポートする
        qiime tools export \
            --input-path ${backbone_fna} \
            --output-path exported_ref

        if [ -f exported_ref/dna-sequences.fasta ]; then
            cp exported_ref/dna-sequences.fasta ${pair_id}/representatives.fasta
        elif ls exported_ref/*.fasta 1> /dev/null 2>&1; then
            cp exported_ref/*.fasta ${pair_id}/representatives.fasta
        else
            cp exported_ref/*.fna ${pair_id}/representatives.fasta
        fi

        # 3. タクソノミ情報をTSVに変換
        qiime tools export \
            --input-path ${taxonomy} \
            --output-path exported_taxonomy

        if [ -f exported_taxonomy/taxonomy.tsv ]; then
            cp exported_taxonomy/taxonomy.tsv ${pair_id}/taxonomy.tsv
        elif [ -f exported_taxonomy/consensus_assignments.tsv ]; then
            cp exported_taxonomy/consensus_assignments.tsv ${pair_id}/taxonomy.tsv
        else
            find exported_taxonomy -name "*.tsv" -exec cp {} ${pair_id}/taxonomy.tsv \\;
        fi

        # 4. feature-table.tsv と taxonomy.tsv を結合し、rrnDBコピー数補正 & 全体和正規化を行い、カウント数順（降順）でソートする (タクソノミ)
        python3 ${params.petagenomeDir}/scripts/Python/parse_taxonomy.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${pair_id}/taxonomy_counts.tsv \
            --mode ${rrndb_mode}

        # 5. アノテーションテーブルを紐づけて、rrnDB補正後の存在量を各機能（KO, MetaCycなど）に分配・集計する
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
workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    input_ch // tuple val(ref_id), val(pair_id), path(sam_file)
    backbone_fna
    taxonomy
    rrndb_stats
    rrndb_mode
    annotation_table
    target_annots

    main:
    out = qiime2_greengenes2_wgs(
        p.combine(input_ch).map { p_val, ref_id, pair_id, sam_file -> tuple(p_val, pair_id, sam_file) },
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
workflow QIIME2_GREENGENES2_WGS_ALL {
    p = createNullParamsChannel()
    mode_ch = Channel.value(params.qiime2_gg2_rrndb_mode)
    
    backbone_ch = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))
    rrndb_ch = Channel.value(file(params.qiime2_gg2_wgs_rrndb_stats, checkIfExists: true))
    annotation_ch = Channel.value(file(params.annotation_table, checkIfExists: true))
    annots_ch = Channel.value(params.functional_annotations)

    wgs_input_ch = Channel.fromPath("${params.output}/MAP_SUB/*/*.sam")
        .map { sam_path ->
            def ref_id = "gg2_backbone"
            def pair_id = sam_path.parent.name
            return tuple(ref_id, pair_id, sam_path)
        }

    QIIME2_GREENGENES2_WGS_SUB(
        p,
        wgs_input_ch,
        backbone_ch,
        taxonomy_ch,
        rrndb_ch,
        mode_ch,
        annotation_ch,
        annots_ch
    )
}

workflow {
    QIIME2_GREENGENES2_WGS_ALL()
}