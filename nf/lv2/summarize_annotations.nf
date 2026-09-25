#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 動的な属性リストのデフォルト値
params.target_attrs = params.containsKey('target_attrs') ? params.target_attrs : ['taxid', 'ko']

// 2. プロセス固有の上限設定 
def SUMMARIZE_MAX_MEMORY  = 8
def SUMMARIZE_MAX_THREADS = 2

// 3. 上限値による動的クリッピング
params.summarize_annotations_memory  = Math.min(params.memory as Integer, SUMMARIZE_MAX_MEMORY)
params.summarize_annotations_threads = Math.min(params.threads as Integer, SUMMARIZE_MAX_THREADS)

// 低相同性アノテーションを除外するためのデフォルト閾値
params.summarize_min_anno_pident = 0.0

include { createNullParamsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. プロセス定義
// ==========================================

process summarize_annotations {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.summarize_annotations_memory}"
    def threads = "${params.summarize_annotations_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(qry_id), path(samples_anno_tsv)
        val target_attrs // 例: ['taxid', 'ko']
        path name_maps, stageAs: 'name_map_*' // 各属性の名前マップ群

    output:
        // 動的に生成されるサマリファイルをタプル（qry_id, [ファイルパスのリスト]）として出力
        tuple val(qry_id), path("${qry_id}.*.summary.tsv"), emit: annotations_summary

    script:
    def py_script  = "${params.petagenomeDir}/scripts/Python/summarize_annotations.py"
    def min_pident = params.summarize_min_anno_pident

    // 各属性に対応するファイルが存在する場合のみ、ステージングされたファイル名とオプションを割り当てる
    def map_opts = ""
    name_maps.eachWithIndex { map_file, idx ->
        if (map_file && map_file.name != 'NO_FILE') {
            def attr = target_attrs[idx]
            map_opts += " --${attr}_name_map ${map_file.name}"
        }
    }

    def attrs_str = target_attrs.join(',')

    """
    echo "${processProfile(task)}" | tee prof.txt

    python3 "${py_script}" \
        -i "${samples_anno_tsv}" \
        -o "${qry_id}" \
        -p "${min_pident}" \
        --target_attrs "${attrs_str}" \
        ${map_opts}
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow SUMMARIZE_ANNOTATIONS_SUB {
    take:
    p
    contig_anno_ch // tuple(qry_id, path_tsv)
    name_maps_map  // Map [attr: path (または null / file('NO_FILE'))]

    main:
    // target_attrs の順序を確定（文字列・リスト両対応）
    def target_attrs = params.target_attrs instanceof Collection ? 
                       params.target_attrs : 
                       params.target_attrs.toString().split(',').collect { it.trim() }

    // あらかじめリストを評価してから Channel.value に包む
    def resolved_maps = target_attrs.collect { attr -> 
        def val = name_maps_map[attr]
        return val ? val : file('NO_FILE')
    }
    def name_maps_ch = Channel.value(resolved_maps)

    // プロセスへ投入
    res = summarize_annotations(contig_anno_ch, target_attrs, name_maps_ch)

    emit:
    annotations_summary = res.annotations_summary
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

workflow SUMMARIZE_ANNOTATIONS_ALL {
    p = createNullParamsChannel()

    // 入力TSVチャンネルの作成
    contig_anno_ch = Channel.fromPath(params.summarize_input_tsv_path)
        .map { f ->
            def qry_id = f.name.replaceAll(/_contig_.*\.tsv$/, '')
            tuple(qry_id, f)
        }

    def target_attrs = params.target_attrs instanceof Collection ? 
                       params.target_attrs : 
                       params.target_attrs.toString().split(',').collect { it.trim() }
    
    def name_maps_map = [:]

    target_attrs.each { attr ->
        def mapPathKey = "${attr}_name_map_path"
        
        // params[...] による直接参照へ変更
        def mapPathVal = params[mapPathKey] ? params[mapPathKey].toString().trim() : null
        name_maps_map[attr] = mapPathVal ? file(mapPathVal, checkIfExists: true) : null
        
        println "DEBUG_CHECK: attr=${attr}, key=${mapPathKey}, val=${mapPathVal}"
    }

    out_ch = SUMMARIZE_ANNOTATIONS_SUB(p, contig_anno_ch, name_maps_map)

    out_ch.annotations_summary.view { qry_id, files -> "ANNOTATIONS SUMMARY TSV for ${qry_id}: $files" }
}

workflow {
    SUMMARIZE_ANNOTATIONS_ALL()
}