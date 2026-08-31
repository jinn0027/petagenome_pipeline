#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 動的な属性リストのデフォルト値
params.target_attrs = params.containsKey('target_attrs') ? params.target_attrs : ['taxid', 'ko']

// 2. プロセス固有の上限設定
def ANNOTATE_SAMPLES_MAX_MEMORY  = 8 // GB (AWK 連想配列でのアノテーション保持用)
def ANNOTATE_SAMPLES_MAX_THREADS = 2 // AWK 処理および I/O の上限

// 3. 上限値による動的クリッピング
params.annotate_samples_memory  = Math.min(params.memory as Integer, ANNOTATE_SAMPLES_MAX_MEMORY)
params.annotate_samples_threads = Math.min(params.threads as Integer, ANNOTATE_SAMPLES_MAX_THREADS)

// 1コンティグ・1ORFあたりの保持上限数のデフォルト値
params.annotate_samples_top_n_samples_vs_merged = 3  // 1 ORF に対して保持する merged ORF の上限
params.annotate_samples_top_n_merged_vs_anno = 3  // 1 merged ORF に対して保持する DB Hit (動的属性) の上限

include { createNullParamsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. プロセス定義
// ==========================================

process annotate_samples {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.annotate_samples_memory}"
    def threads = "${params.annotate_samples_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(qry_id), path(samples_map_m8), path(annotated_tsv)
        val target_attrs                     // 例: ['taxid', 'ko'] または 'taxid,ko'

    output:
        tuple val(qry_id), path("${qry_id}_samples_annotated.tsv"), emit: samples_anno

    script:
    def map_top_n  = params.annotate_samples_top_n_samples_vs_merged
    def anno_top_n = params.annotate_samples_top_n_merged_vs_anno
    
    // 入力が文字列の場合も考慮して確実にリスト化してから空白区切り文字列にする
    def attrs_list = target_attrs instanceof Collection ? 
                       target_attrs : 
                       target_attrs.toString().split(',').collect { it.trim() }
    def attrs_str  = attrs_list.join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Bash用のシェル変数として明示的に保持（unbound variable エラー防止）
    TARGET_ATTRS_STR="${attrs_str}"

    awk -v ANNO_TOP_N="${anno_top_n}" \
        -v MAP_TOP_N="${map_top_n}" \
        -v TARGET_ATTRS="\${TARGET_ATTRS_STR}" '
    BEGIN {
        FS="\\t"; OFS="\\t"
        n_attrs = split(TARGET_ATTRS, attr_list, " ")
    }

    # -------------------------------------------------------------
    # 1. アノテーション結果 (annotated.tsv) の読み込み
    # -------------------------------------------------------------
    NR == FNR {
        if (FNR == 1) {
            for (i = 1; i <= NF; i++) {
                col[\$i] = i
            }
            c_qseqid   = ("qseqid" in col   ? col["qseqid"]   : 1)
            c_sseqid   = ("sseqid" in col   ? col["sseqid"]   : 2)
            c_pident   = ("pident" in col   ? col["pident"]   : 3)
            c_evalue   = ("evalue" in col   ? col["evalue"]   : 11)
            c_bitscore = ("bitscore" in col ? col["bitscore"] : 12)
            
            for (k = 1; k <= n_attrs; k++) {
                an = attr_list[k]
                c_attrs[k] = (an in col ? col[an] : (12 + k))
            }
            next
        }

        orf_id   = \$c_qseqid
        sseqid   = \$c_sseqid
        pident   = \$c_pident
        evalue   = \$c_evalue
        bitscore = \$c_bitscore
        
        all_na = 1
        for (k = 1; k <= n_attrs; k++) {
            val = (NF >= c_attrs[k] && \$c_attrs[k] != "" && \$c_attrs[k] != " ") ? \$c_attrs[k] : "N/A"
            current_attr_val[k] = val
            if (val != "N/A" && val != "") {
                all_na = 0
            }
        }

        if (all_na) {
            next
        }

        if (anno_count[orf_id] < ANNO_TOP_N) {
            anno_count[orf_id]++
            idx = anno_count[orf_id]

            orf_sseqid[orf_id, idx]   = (sseqid   != "" ? sseqid   : "N/A")
            orf_pident[orf_id, idx]   = (pident   != "" ? pident   : "N/A")
            orf_evalue[orf_id, idx]   = (evalue   != "" ? evalue   : "N/A")
            orf_bitscore[orf_id, idx] = (bitscore != "" ? bitscore : "N/A")
            
            for (k = 1; k <= n_attrs; k++) {
                orf_attr_vals[orf_id, idx, k] = current_attr_val[k]
            }
        }
        next
    }

    # -------------------------------------------------------------
    # 2. samples -> ORF マッピング結果 (m8) の処理
    # -------------------------------------------------------------
    {
        orf_in_sample_id = \$1
        orf_id     = \$2
        c_pident   = \$3
        c_evalue   = \$11
        c_bitscore = \$12

        if (map_count[orf_in_sample_id] < MAP_TOP_N) {
            map_count[orf_in_sample_id]++

            if (orf_id in anno_count) {
                n = anno_count[orf_id]
                for (i = 1; i <= n; i++) {
                    # 出力行の構築
                    out_str = orf_in_sample_id FS orf_id FS c_pident FS c_evalue FS c_bitscore FS orf_sseqid[orf_id, i]
                    for (k = 1; k <= n_attrs; k++) {
                        out_str = out_str FS orf_attr_vals[orf_id, i, k]
                    }
                    out_str = out_str FS orf_pident[orf_id, i] FS orf_evalue[orf_id, i] FS orf_bitscore[orf_id, i]
                    print out_str
                }
            } else {
                out_str = orf_in_sample_id FS orf_id FS c_pident FS c_evalue FS c_bitscore FS "N/A"
                for (k = 1; k <= n_attrs; k++) {
                    out_str = out_str FS "N/A"
                }
                out_str = out_str FS "N/A" FS "N/A" FS "N/A"
                print out_str
            }
        }
    }
    ' "${annotated_tsv}" "${samples_map_m8}" > temp_body.tsv

    # 動的ヘッダーを生成してファイルに書き出す (TARGET_ATTRS_STR を使用)
    HEADER="orf_in_sample_id\\torf_id\\tmap_pident\\tmap_evalue\\tmap_bitscore\\tref_sseqid"
    for attr in \${TARGET_ATTRS_STR}; do
        HEADER="\${HEADER}\\t\${attr}"
    done
    HEADER="\${HEADER}\\tanno_pident\\tanno_evalue\\tanno_bitscore"

    echo -e "\${HEADER}" > "${qry_id}_samples_annotated.tsv"
    cat temp_body.tsv >> "${qry_id}_samples_annotated.tsv"
    rm -f temp_body.tsv
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow ANNOTATE_SAMPLES_SUB {
    take:
    p
    samples_map_out // tuple(qry_id, path_m8)
    annotated_res   // path(tsv) または tuple(ref_id, path_tsv)

    main:
    samples_map_out.view {
        "[DEBUG VIEW] ANNOTATE_SAMPLES_SUB: samples_map_out\n" +
        "  - 期待される形式: tuple(val(qry_id), path(m8_file))\n" +
        "  - 実際の値     : ${it}"
    }

    annotated_res.view {
        "[DEBUG VIEW] ANNOTATE_SAMPLES_SUB: annotated_res\n" +
        "  - 期待される形式: path(tsv_file) または tuple(val(id), path(tsv_file))\n" +
        "  - 実際の値     : ${it}"
    }

    annotated_path_ch = annotated_res.map { it instanceof Path ? it : it[1] }

    joined_ch = samples_map_out
        .combine(annotated_path_ch.collect())
        .map { qry_id, map_m8, anno_tsv ->
            tuple(qry_id, map_m8, anno_tsv)
        }

    joined_ch.view {
        "[DEBUG VIEW] ANNOTATE_SAMPLES_SUB: joined_ch (Process Input)\n" +
        "  - 実際の値     : ${it}"
    }

    // カンマ区切りの文字列で渡された場合もリストとしてパース
    def target_attrs = params.target_attrs instanceof Collection ? 
                       params.target_attrs : 
                       params.target_attrs.toString().split(',').collect { it.trim() }
    
    res = annotate_samples(joined_ch, target_attrs)

    res.samples_anno.view {
        "[DEBUG VIEW] ANNOTATE_SAMPLES_SUB: res.samples_anno (Output)\n" +
        "  - 実際の値     : ${it}"
    }

    emit:
    samples_anno = res.samples_anno
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

workflow ANNOTATE_SAMPLES_ALL {
    p = createNullParamsChannel()

    samples_map_ch = Channel.fromPath(params.annotate_samples_m8_path)
                          .map { f -> 
                              def qry_id = f.name.replaceAll(/_samples_map\.m8$/, '')
                              tuple(qry_id, f) 
                          }
    annotated_ch  = Channel.fromPath(params.annotate_samples_tsv_path)
                          .map { f -> tuple('merged_all_samples', f) }

    out_ch = ANNOTATE_SAMPLES_SUB(p, samples_map_ch, annotated_ch)
    out_ch.samples_anno.view { i -> "SAMPLES ANNO TSV: $i" }
}

workflow {
    ANNOTATE_SAMPLES_ALL()
}