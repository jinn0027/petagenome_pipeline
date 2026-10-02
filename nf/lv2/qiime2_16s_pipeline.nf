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
params.qiime2_gg2_perc_identity = 0.99 // パーセント一致度閾値のデフォルト設定
params.qiime2_gg2_rrndb_mode = 'right' // rrnDBコピー数探索モードのデフォルト
params.qiime2_gg2_only_gxxxx = true    // G[0-9]ノードのみに絞る（デフォルト: true）

// 機能アノテーション関連パラメータのデフォルト追加
if (!params.containsKey('annotation_table')) {
    params.annotation_table = "${params.petagenomeDir}/data/greengenes2/genome_annotations_table.tsv"
}
params.functional_annotations = "KO,MetaCyc"

params.qiime2_gg2_16s_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_16s_taxonomy = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_gg2_16s_rrndb_stats = "${params.petagenomeDir}/data/rrnDB/rrnDB-5.10_pantaxa_stats_RDP.tsv.gz"

// 必須パラメータのチェック
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it."
}

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; createSeqsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// 下位モジュール（lv1）のインポート
include { QIIME2_DADA2_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_dada2"
include { QIIME2_GREENGENES2_16S_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_gg2_16s"

// ==========================================
// バックボーンを G[0-9] のみにフィルタリングするプロセス (16S用)
// ==========================================
process filter_backbone_gxxxx_16s {
    tag "filtering backbone for G[0-9] (16S)"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }

    def gb = "${params.qiime2_gg2_16s_memory ?: params.memory}"
    def threads = "${params.qiime2_gg2_16s_threads ?: params.threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
    path backbone_fna

    output:
    path "filtered_backbone.fna.qza", emit: backbone

    script:
    """
    export PYTHONWARNINGS="ignore"
    export XDG_CONFIG_HOME=/tmp/qiime2_config
    export MPLCONFIGDIR=/tmp/matplotlib_config
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    export FONTCONFIG_PATH=/tmp/fontconfig

    qiime tools export \
        --input-path ${backbone_fna} \
        --output-path exported_backbone

    python3 - << 'EOF'
import re

input_file = "exported_backbone/dna-sequences.fasta"
output_file = "filtered_backbone.fasta"
pattern = re.compile(r"^G[0-9]")

count = 0
filtered_count = 0

with open(input_file, "r") as fin, open(output_file, "w") as fout:
    write_this = False
    for line in fin:
        if line.startswith(">"):
            count += 1
            header = line.strip()
            seq_id = header[1:].split()[0]
            if pattern.match(seq_id):
                write_this = True
                filtered_count += 1
                print(header, file=fout)
            else:
                write_this = False
        else:
            if write_this:
                fout.write(line)

print(f"Total backbone sequences: {count}")
print(f"Filtered G[0-9] sequences: {filtered_count}")
EOF

    qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path filtered_backbone.fasta \
        --output-path filtered_backbone.fna.qza
    """
}

// ==========================================
// 追加：FASTA直接入力からDADA2出力を模倣するプロセス（標準ライブラリのみ使用）
// ==========================================
process qiime2_import_fasta_as_features {
    tag "${sample_id}"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    publishDir "${params.output}/qiime2_dada2/${sample_id}", mode: 'symlink', enabled: params.publish_output

    input:
        tuple val(sample_id), path(fasta_file)

    output:
        tuple val(sample_id), 
              path("${sample_id}/table.qza"), 
              path("${sample_id}/rep-seqs.qza"), 
              path("${sample_id}/denoising-stats.qza"), 
              path("${sample_id}/base-transition-stats.qza")

    script:
        """
        export PYTHONWARNINGS="ignore"
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config

        mkdir -p ${sample_id}

        # 1. 代表配列 (rep-seqs.qza) のインポート
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path ${fasta_file} \
            --output-path ${sample_id}/rep-seqs.qza

        # 2. 標準PythonだけでFASTAをパースしてBIOMテーブルを生成
        python3 - <<EOF
import biom
import pandas as pd

fasta_path = "${fasta_file}"
record_ids = []

with open(fasta_path, 'r') as f:
    for line in f:
        if line.startswith('>'):
            rec_id = line[1:].strip().split()[0]
            record_ids.append(rec_id)

counts = list(range(1, len(record_ids) + 1))
data = pd.DataFrame(counts, index=record_ids, columns=["${sample_id}"])
table = biom.Table(data.values, data.index, data.columns)

biom_file = "${sample_id}/feature_table.biom"
with open(biom_file, "w") as f:
    table.to_json("Generated by pipeline fasta bypass", f)
EOF

        # 3. フィーチャーテーブル (table.qza) のインポート
        qiime tools import \
            --type 'FeatureTable[Frequency]' \
            --input-path ${sample_id}/feature_table.biom \
            --input-format BIOMV100Format \
            --output-path ${sample_id}/table.qza

        touch ${sample_id}/denoising-stats.qza
        touch ${sample_id}/base-transition-stats.qza
        """
}

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
    rrndb_stats
    rrndb_mode
    annotation_table
    target_annots
    perc_identity
    is_fasta

    main:
    // FASTAの場合はDADA2をスキップし、それ以外は従来のDADA2をそのまま実行
    if (is_fasta) {
        dada2_out = qiime2_import_fasta_as_features(reads)
    } else {
        dada2_out = QIIME2_DADA2_SUB(p, reads)
    }

    // B-0. パラメータに応じてバックボーンを G[0-9] ノードのみにフィルタリング
    ch_backbone = params.qiime2_gg2_only_gxxxx ? filter_backbone_gxxxx_16s(backbone_fna) : backbone_fna

    // B. Greengenes2 (16S用) による系統配置・タクソノミ付与 & 機能アノテーション集計
    gg2_out = QIIME2_GREENGENES2_16S_SUB(
        p,
        dada2_out,
        target_region,
        ch_backbone, // フィルタ済み（または元の）バックボーンを渡す
        taxonomy,
        rrndb_stats,
        rrndb_mode,
        annotation_table,
        target_annots,
        perc_identity
    )

    emit:
    gg2_out = gg2_out
}

// ==========================================
// コマンドライン用エントリーポイント (-entry)
// ==========================================
workflow QIIME2_16S_PIPELINE_ALL {
    p = createNullParamsChannel()

    def is_fasta = false
    if (params.containsKey('qiime2_reads')) {
        def r_str = params.qiime2_reads.toString().toLowerCase()
        if (r_str.endsWith('.fasta') || r_str.endsWith('.fa') || r_str.endsWith('.fna')) {
            is_fasta = true
        }
    }

    def reads
    if (is_fasta) {
        reads = createSeqsChannel(params.qiime2_reads)
    } else {
        reads = createPairsChannel(params.qiime2_reads)
    }

    region_ch = Channel.value(params.qiime2_gg2_16s_target_region)
    mode_ch = Channel.value(params.qiime2_gg2_rrndb_mode)
    perc_identity_ch = Channel.value(params.qiime2_gg2_perc_identity)
    
    annotation_ch = Channel.value(file(params.annotation_table, checkIfExists: true))
    annots_ch = Channel.value(params.functional_annotations)

    backbone_ch = Channel.value(file(params.qiime2_gg2_16s_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_16s_taxonomy, checkIfExists: true))
    rrndb_ch = Channel.value(file(params.qiime2_gg2_16s_rrndb_stats, checkIfExists: true))

    out_ch = QIIME2_16S_PIPELINE_SUB(
        p,
        reads,
        region_ch,
        backbone_ch,
        taxonomy_ch,
        rrndb_ch,
        mode_ch,
        annotation_ch,
        annots_ch,
        perc_identity_ch,
        is_fasta
    )

    out_ch.gg2_out.view { i -> "QIIME2 16S PIPELINE OUT: $i" }
}

workflow {
    QIIME2_16S_PIPELINE_ALL()
}