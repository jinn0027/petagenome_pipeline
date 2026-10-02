#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_WGS_MAX_MEMORY = 64
def QIIME2_WGS_MAX_THREADS = 32

params.qiime2_wgs_memory = Math.min(params.memory as Integer, QIIME2_WGS_MAX_MEMORY)
params.qiime2_wgs_threads = Math.min(params.threads as Integer, QIIME2_WGS_MAX_THREADS)

// パラメータの初期化
params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy     = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_gg2_wgs_rrndb_stats  = "${params.petagenomeDir}/data/rrnDB/rrnDB-5.10_pantaxa_stats_RDP.tsv.gz"
params.qiime2_gg2_rrndb_mode       = 'right' // rrnDBコピー数探索モードのデフォルト
params.qiime2_gg2_only_gxxxx       = true    // G[0-9]ノードのみに絞る（デフォルト: true）

// 機能アノテーション関連パラメータのデフォルト追加
if (!params.containsKey('annotation_table')) {
    params.annotation_table = "${params.petagenomeDir}/data/greengenes2/genome_annotations_table.tsv"
}
params.functional_annotations = "KO,MetaCyc"

// 必須パラメータのチェック
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it."
}

include { createNullParamsChannel; createPairsChannel; apptainerContainerOptions } from "${params.petagenomeDir}/nf/common/utils"

// 下位モジュールのインポート
include { FASTP_SUB } from "${params.petagenomeDir}/nf/lv1/fastp"
include { BUILD_REF_DB_SUB; MAP_SUB } from "${params.petagenomeDir}/nf/lv1/bowtie2"
include { QIIME2_GREENGENES2_WGS_SUB } from "${params.petagenomeDir}/nf/lv1/qiime2_gg2_wgs"

// ==========================================
// バックボーンを G[0-9] のみにフィルタリングするプロセス (WGS用)
// ==========================================
process filter_backbone_gxxxx_wgs {
    tag "filtering backbone for G[0-9] (WGS)"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }

    def gb = "${params.qiime2_wgs_memory ?: params.memory}"
    def threads = "${params.qiime2_wgs_threads ?: params.threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads

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
    rrndb_stats
    rrndb_mode
    annotation_table
    target_annots

    main:

    // A. fastp による品質管理・フィルタリング
    fastp_out = FASTP_SUB(p, reads)

    // B-0. パラメータに応じてバックボーンを G[0-9] ノードのみにフィルタリング
    ch_backbone = params.qiime2_gg2_only_gxxxx ? filter_backbone_gxxxx_wgs(backbone) : backbone

    // B-1. フィルタ済み .qza をプレーンな FASTA にエクスポート
    exported_ref = export_qiime2_fasta(ch_backbone)
    exported_ref_ch = exported_ref.map { fasta_path ->
        def ref_id = "gg2_backbone"
        return tuple(ref_id, [fasta_path])
    }
    
    // B-2. Bowtie2 リファレンス DB の作成
    ref_db = BUILD_REF_DB_SUB(p, exported_ref_ch)

    // B-3. Bowtie2 によるリファレンスへのマッピング
    bowtie2_out = MAP_SUB(p, ref_db, fastp_out)

    // C. Greengenes2 によるタクソノミ付与・集計 & 機能アノテーション集計
    gg2_out = QIIME2_GREENGENES2_WGS_SUB(
        p, 
        bowtie2_out, 
        ch_backbone, 
        taxonomy,
        rrndb_stats,
        rrndb_mode,
        annotation_table,
        target_annots
    )

    emit:
    gg2_out = gg2_out
}

// ==========================================
// コマンドライン用エントリーポイント (-entry)
// ==========================================
workflow QIIME2_WGS_PIPELINE_ALL {
    p         = createNullParamsChannel()
    reads     = createPairsChannel(params.qiime2_reads)
    mode_ch   = Channel.value(params.qiime2_gg2_rrndb_mode)
    
    // 機能アノテーション関連チャンネルの作成
    annotation_ch = Channel.value(file(params.annotation_table, checkIfExists: true))
    annots_ch     = Channel.value(params.functional_annotations)

    // 必要なリファレンスファイルのロード
    backbone_ch  = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch  = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))
    rrndb_ch     = Channel.value(file(params.qiime2_gg2_wgs_rrndb_stats, checkIfExists: true))

    out_ch = QIIME2_WGS_PIPELINE_SUB(
        p,
        reads,
        backbone_ch,
        taxonomy_ch,
        rrndb_ch,
        mode_ch,
        annotation_ch,
        annots_ch
    )

    out_ch.gg2_out.view { i -> "QIIME2 WGS PIPELINE OUT: $i" }
}

workflow {
    QIIME2_WGS_PIPELINE_ALL()
}