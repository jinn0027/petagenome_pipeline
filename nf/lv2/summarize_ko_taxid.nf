#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 低相同性アノテーションを除外するためのデフォルト閾値 (例: 0.0)
params.summarize_min_anno_pident = 0.0

include { createNullParamsChannel; clusterOptions; processProfile } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. プロセス定義
// ==========================================

process SUMMARIZE_KO_TAXID {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    // Pythonスクリプトが配置されているディレクトリを参照できるよう petagenomeDir をマウント(bind)
    containerOptions = "${params.apptainerRunOptions} -B ${params.petagenomeDir}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.memory ?: 4}"
    def threads = "${params.cpus ?: 1}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
    tuple val(qry_id), path(contig_anno_tsv)
    path ko_name_map, stageAs: 'ko_map.tsv'     // params.ko_name_map から渡す (任意)
    path taxid_name_map, stageAs: 'taxid_map.tsv'  // params.taxid_name_map から渡す (任意)

    output:
    tuple val(qry_id), path("${qry_id}.ko_summary.tsv")   , emit: ko_summary
    tuple val(qry_id), path("${qry_id}.taxid_summary.tsv"), emit: tax_summary

    script:
    def py_script   = "${params.petagenomeDir}/scripts/Python/summarize_ko_taxid.py"
    def min_pident  = params.summarize_min_anno_pident
    
    // stageAs で指定されたファイル名を用いてオプションを組み立て
    def ko_opt  = (ko_name_map.name != 'NO_FILE') ? "--ko_name_map ko_map.tsv" : ""
    def tax_opt = (taxid_name_map.name != 'NO_FILE') ? "--taxid_name_map taxid_map.tsv" : ""

    """
    #!/usr/bin/env bash
    set -euo pipefail

    python3 "${py_script}" \\
        -i "${contig_anno_tsv}" \\
        -o "${qry_id}" \\
        -p "${min_pident}" \\
        ${ko_opt} \\
        ${tax_opt}
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow SUMMARIZE_KO_TAXID_SUB {
    take:
    p
    contig_anno_ch // tuple(qry_id, path_tsv)
    ko_name_map    // path (または file('NO_FILE'))
    taxid_name_map // path (または file('NO_FILE'))

    main:
    res = SUMMARIZE_KO_TAXID(contig_anno_ch, ko_name_map, taxid_name_map)

    emit:
    ko_summary  = res.ko_summary
    tax_summary = res.tax_summary
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

workflow SUMMARIZE_KO_TAXID_ALL {
    p = createNullParamsChannel()

    // 入力TSVチャンネルの作成 (*_contig_taxid_ko.tsv を想定)
    contig_anno_ch = Channel.fromPath(params.summarize_input_tsv_path)
        .map { f ->
            def qry_id = f.name.replaceAll(/_contig_taxid_ko\.tsv\$/, '')
            tuple(qry_id, f)
        }

    // パラメータ指定の有無を判定してファイルパスを取得
    ko_map_file  = (params.containsKey('ko_name_map_path') && params.ko_name_map_path) ? file(params.ko_name_map_path, checkIfExists: true) : file('NO_FILE')
    tax_map_file = (params.containsKey('taxid_name_map_path') && params.taxid_name_map_path) ? file(params.taxid_name_map_path, checkIfExists: true) : file('NO_FILE')

    out_ch = SUMMARIZE_KO_TAXID_SUB(p, contig_anno_ch, ko_map_file, tax_map_file)
    
    out_ch.ko_summary.view { i -> "KO SUMMARY TSV: $i" }
    out_ch.tax_summary.view { i -> "TAXID SUMMARY TSV: $i" }
}

workflow {
    SUMMARIZE_KO_TAXID_ALL()
}