#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_WOLTKA_MAX_MEMORY  = 64
def QIIME2_WOLTKA_MAX_THREADS = 32

params.qiime2_woltka_memory  = Math.min(params.memory as Integer, QIIME2_WOLTKA_MAX_MEMORY)
params.qiime2_woltka_threads = Math.min(params.threads as Integer, QIIME2_WOLTKA_MAX_THREADS)

// リファレンスファイルのパス（未定義・未指定の場合は即座にエラーにする）
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_woltka_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_woltka_taxonomy     = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_woltka_wol_map      = "${params.petagenomeDir}/data/greengenes2/taxonomy/curated/taxid/taxid.map"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_woltka {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_woltka_memory}"
    def threads = "${params.qiime2_woltka_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(ref_id), val(qry_id), path(sam), path(wol_map), path(backbone), path(taxonomy)

    output:
        tuple val(ref_id), val(qry_id), 
              path("${qry_id}/feature-table.tsv"), 
              path("${qry_id}/taxonomy.tsv")

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
        mkdir -p ${qry_id}

        # 1. Woltkaによるアライメント結果（SAM）からのプロファイル集計
        woltka classify \\
            --input ${sam} \\
            --map ${wol_map} \\
            --output ${qry_id}/woltka_table.biom \\
            --threads ${threads}

        # 2. 生成されたBIOMテーブルをQIIME 2アーティファクト（.qza）にインポート
        qiime tools import \\
            --type 'FeatureTable[Frequency]' \\
            --input-path ${qry_id}/woltka_table.biom \\
            --output-path ${qry_id}/table.qza

        # 3. Greengenes2 バックボーンへのマッピング・同期
        qiime greengenes2 non-v4-16s \\
            --i-table ${qry_id}/table.qza \\
            --i-backbone ${backbone} \\
            --o-mapped-table ${qry_id}/mapped_table.qza \\
            --o-representatives ${qry_id}/representatives.qza \\
            --p-threads ${threads}

        # 4. フィーチャーテーブルをTSV（テキスト）に変換
        qiime tools export \\
            --input-path ${qry_id}/mapped_table.qza \\
            --output-path ${qry_id}/exported_table
        
        biom convert \\
            -i ${qry_id}/exported_table/feature-table.biom \\
            -o ${qry_id}/feature-table.tsv \\
            --to-tsv

        # 5. タクソノミー情報をTSVに変換
        qiime tools export \\
            --input-path ${taxonomy} \\
            --output-path ${qry_id}/exported_taxonomy

        if [ -f ${qry_id}/exported_taxonomy/taxonomy.tsv ]; then
            cp ${qry_id}/exported_taxonomy/taxonomy.tsv ${qry_id}/taxonomy.tsv
        elif [ -f ${qry_id}/exported_taxonomy/consensus_assignments.tsv ]; then
            cp ${qry_id}/exported_taxonomy/consensus_assignments.tsv ${qry_id}/taxonomy.tsv
        else
            for f in ${qry_id}/exported_taxonomy/*.tsv; do
                if [ -f "\$f" ]; then
                    cp "\$f" ${qry_id}/taxonomy.tsv
                    break
                fi
            done
        fi
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow QIIME2_WOLTKA_SUB {
    take:
    p
    map_out        
    wol_map        
    backbone       
    taxonomy       

    main:
    // map_out をベースにして結合する際、バリューチャンネルの「値(.val)」を確実に抽出して渡す
    in_ch = map_out.combine(p)
        .map { ref_id, qry_id, sam, p_val ->
            // バリューチャンネルからファイル／値オブジェクトを .val で安全に取り出す
            return tuple(
                p_val, 
                ref_id, 
                qry_id, 
                sam, 
                wol_map.val,   // ← .val を付与して中身のファイルを取り出す
                backbone.val,  // ← 同上
                taxonomy.val   // ← 同上
            )
        }

    out = qiime2_woltka(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン用エントリーポイント
// ==========================================
workflow QIIME2_WOLTKA_ALL {
    p            = createNullParamsChannel()
    
    wol_map_ch   = Channel.value(file(params.qiime2_woltka_wol_map, checkIfExists: true))
    backbone_ch  = Channel.value(file(params.qiime2_woltka_backbone_fna, checkIfExists: true))
    taxonomy_ch  = Channel.value(file(params.qiime2_woltka_taxonomy, checkIfExists: true))

    map_dummy_ch = Channel.fromPath("${params.output}/bowtie2/*/out.sam")
        .map { sam_path ->
            def qry_id = sam_path.parent.name
            def ref_id = sam_path.parent.parent.name
            return tuple(ref_id, qry_id, sam_path)
        }

    QIIME2_WOLTKA_SUB(
        p,
        map_dummy_ch,
        wol_map_ch,
        backbone_ch,
        taxonomy_ch
    )
}

workflow {
    QIIME2_WOLTKA_ALL()
}
