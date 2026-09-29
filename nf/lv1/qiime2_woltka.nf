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

if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_woltka_wol_map     = "${params.petagenomeDir}/data/greengenes2/taxonomy/curated/taxid/taxid.map"
params.qiime2_woltka_taxonomy    = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_woltka_target_rank = "species"

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
        tuple val(p), val(ref_id), val(qry_id), path(sam), path(wol_map), path(taxonomy)

    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/table.qza")

    script:
        def targetRank = params.qiime2_woltka_target_rank
        """
        export PYTHONWARNINGS="ignore"
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${qry_id}

        qiime tools import \
            --type 'FeatureData[SeqAlnMap]' \
            --input-path ${sam} \
            --output-path ${qry_id}/alignment.qza

        qiime tools import \
            --type 'FeatureData[SimpleMap]' \
            --input-path ${wol_map} \
            --output-path ${qry_id}/taxon_map.qza

        qiime woltka classify \
            --i-alignment ${qry_id}/alignment.qza \
            --i-taxon-map ${qry_id}/taxon_map.qza \
            --i-reference-taxonomy ${taxonomy} \
            --p-target-rank ${targetRank} \
            --o-classified-table ${qry_id}/table.qza
        """
}

workflow QIIME2_WOLTKA_SUB {
    take:
    p
    map_out        
    wol_map        
    taxonomy       

    main:
    in_ch = map_out.combine(p)
        .map { ref_id, qry_id, sam, p_val ->
            return tuple(p_val, ref_id, qry_id, sam, wol_map.val, taxonomy.val)
        }

    out = qiime2_woltka(in_ch)
    emit:
    out = out
}

workflow QIIME2_WOLTKA_ALL {
    p            = createNullParamsChannel()
    wol_map_ch   = Channel.value(file(params.qiime2_woltka_wol_map, checkIfExists: true))
    taxonomy_ch  = Channel.value(file(params.qiime2_woltka_taxonomy, checkIfExists: true))

    map_dummy_ch = Channel.fromPath("${params.output}/bowtie2/*/out.sam")
        .map { sam_path ->
            def qry_id = sam_path.parent.name
            def ref_id = sam_path.parent.parent.name
            return tuple(ref_id, qry_id, sam_path)
        }

    QIIME2_WOLTKA_SUB(p, map_dummy_ch, wol_map_ch, taxonomy_ch)
}

workflow {
    QIIME2_WOLTKA_ALL()
}