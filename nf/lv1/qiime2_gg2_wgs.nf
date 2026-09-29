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

if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

// 系統樹ファイル（.nwk）を指定
params.qiime2_gg2_wgs_tree = "${params.petagenomeDir}/data/greengenes2/2024.09.taxonomy.id.nwk"

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
        path tree_file

    output:
        tuple val(pair_id), 
              path("${pair_id}/feature-table.tsv"), 
              path("${pair_id}/taxonomy.tsv")

    script:
        """
        export PYTHONWARNINGS="ignore"
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # 1. 系統樹ファイル（.nwk）を QIIME 2 アーティファクトとしてインポート
        qiime tools import \
            --type 'Phylogeny[Rooted]' \
            --input-path ${tree_file} \
            --output-path reference_tree.qza

        # 2. 系統樹から、テーブルのフィーチャーに対応する Greengenes2 体系のタクソノミー情報を引き出す
        qiime greengenes2 taxonomy-from-table \
            --i-table ${table} \
            --i-reference-taxonomy reference_tree.qza \
            --o-classification taxonomy.qza

        # 3. フィーチャーテーブルをTSVに変換
        qiime tools export \
            --input-path ${table} \
            --output-path exported_table
        
        biom convert \
            -i exported_table/feature-table.biom \
            -o ${pair_id}/feature-table.tsv \
            --to-tsv

        # 4. Greengenes2 体系に紐付けられたタクソノミー情報をTSVに変換
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

workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    woltka_out   
    tree_file    

    main:
    in_ch = woltka_out.map { pair_id, table ->
        tuple(pair_id, table)
    }

    out = qiime2_greengenes2_wgs(
        p.combine(in_ch).map { p_val, pair_id, table -> tuple(p_val, pair_id, table) },
        tree_file
    )

    emit:
    out = out
}

workflow QIIME2_GREENGENES2_WGS_ALL {
    p           = createNullParamsChannel()
    tree_ch     = Channel.value(file(params.qiime2_gg2_wgs_tree, checkIfExists: true))

    woltka_dummy_ch = Channel.fromPath("${params.output}/qiime2_woltka/*/table.qza")
        .map { table_path ->
            def pair_id = table_path.parent.name
            return tuple(pair_id, table_path)
        }

    QIIME2_GREENGENES2_WGS_SUB(
        p,
        woltka_dummy_ch,
        tree_ch
    )
}

workflow {
    QIIME2_GREENGENES2_WGS_ALL()
}