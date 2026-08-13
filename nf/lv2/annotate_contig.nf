#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def ANNOTATE_CONTIG_MAX_MEMORY  = 8 // GB (AWK 連想配列でのアノテーション保持用)
def ANNOTATE_CONTIG_MAX_THREADS = 2 // AWK 処理および I/O の上限

// 3. 上限値による動的クリッピング
params.annotate_contig_memory  = Math.min(params.memory as Integer, ANNOTATE_CONTIG_MAX_MEMORY)
params.annotate_contig_threads = Math.min(params.threads as Integer, ANNOTATE_CONTIG_MAX_THREADS)

// 1コンティグ・1ORFあたりの保持上限数のデフォルト値
params.annotate_contig_top_n_contig = 3  // 1 ORF に対して保持する merged ORF の上限
params.annotate_contig_top_n_anno   = 3  // 1 merged ORF に対して保持する DB Hit (TaxID/KO) の上限

include { createNullParamsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. プロセス定義
// ==========================================

process assign_taxid_ko_to_contigs {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.annotate_contig_memory}"
    def threads = "${params.annotate_contig_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        // 受け取るタプルを [ qry_id, contig_map_m8, annotated_tsv ] の 3 要素に正確に適合
        tuple val(qry_id), path(contig_map_m8), path(annotated_tsv)

    output:
        tuple val(qry_id), path("${qry_id}_contig_taxid_ko.tsv"), emit: contig_anno

    script:
    def map_top_n  = params.annotate_contig_top_n_contig
    def anno_top_n = params.annotate_contig_top_n_anno

    """
    #!/usr/bin/env bash
    set -euo pipefail

    awk -v ANNO_TOP_N="${anno_top_n}" -v MAP_TOP_N="${map_top_n}" '
    BEGIN { FS="\\t"; OFS="\\t" }

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
            c_taxid    = ("taxid" in col    ? col["taxid"]    : 13)
            c_ko       = ("ko" in col       ? col["ko"]       : 14)
            next
        }

        orf_id   = \$c_qseqid
        sseqid   = \$c_sseqid
        pident   = \$c_pident
        evalue   = \$c_evalue
        bitscore = \$c_bitscore
        
        taxid    = (NF >= c_taxid  && \$c_taxid  != "" && \$c_taxid  != " ") ? \$c_taxid  : "N/A"
        ko       = (NF >= c_ko     && \$c_ko     != "" && \$c_ko     != " ") ? \$c_ko     : "N/A"

        if ((taxid == "N/A" || taxid == "") && (ko == "N/A" || ko == "")) {
            next
        }

        if (anno_count[orf_id] < ANNO_TOP_N) {
            anno_count[orf_id]++
            idx = anno_count[orf_id]

            orf_sseqid[orf_id, idx]   = (sseqid   != "" ? sseqid   : "N/A")
            orf_pident[orf_id, idx]   = (pident   != "" ? pident   : "N/A")
            orf_evalue[orf_id, idx]   = (evalue   != "" ? evalue   : "N/A")
            orf_bitscore[orf_id, idx] = (bitscore != "" ? bitscore : "N/A")
            orf_taxid[orf_id, idx]    = taxid
            orf_ko[orf_id, idx]       = ko
        }
        next
    }

    # -------------------------------------------------------------
    # 2. Contig -> ORF マッピング結果 (m8) の処理
    # -------------------------------------------------------------
    {
        contig_id  = \$1
        orf_id     = \$2
        c_pident   = \$3
        c_evalue   = \$11
        c_bitscore = \$12

        if (map_count[contig_id] < MAP_TOP_N) {
            map_count[contig_id]++

            if (orf_id in anno_count) {
                n = anno_count[orf_id]
                for (i = 1; i <= n; i++) {
                    print contig_id, \
                          orf_id, \
                          c_pident, \
                          c_evalue, \
                          c_bitscore, \
                          orf_sseqid[orf_id, i], \
                          orf_taxid[orf_id, i], \
                          orf_ko[orf_id, i], \
                          orf_pident[orf_id, i], \
                          orf_evalue[orf_id, i], \
                          orf_bitscore[orf_id, i]
                }
            } else {
                print contig_id, \
                      orf_id, \
                      c_pident, \
                      c_evalue, \
                      c_bitscore, \
                      "N/A", \
                      "N/A", \
                      "N/A", \
                      "N/A", \
                      "N/A", \
                      "N/A"
            }
        }
    }
    ' "${annotated_tsv}" "${contig_map_m8}" > temp_body.tsv

    # ヘッダーを付与して出力
    echo -e "contig_id\\torf_id\\tmap_pident\\tmap_evalue\\tmap_bitscore\\tref_sseqid\\ttaxid\\tko\\tanno_pident\\tanno_evalue\\tanno_bitscore" > "${qry_id}_contig_taxid_ko.tsv"
    cat temp_body.tsv >> "${qry_id}_contig_taxid_ko.tsv"
    rm -f temp_body.tsv
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow ANNOTATE_CONTIG_SUB {
    take:
    p
    contig_map_out // tuple(qry_id, path_m8)        <- 2要素
    annotated_res  // tuple(ref_id, path_tsv)      <- 2要素

    main:
    // 2要素 + 2要素 の combine（計4要素）を受け取り、
    // プロセスへ渡す 3要素 [ qry_id, map_m8, anno_tsv ] へ正確にマッピング
    joined_ch = contig_map_out
        .combine(annotated_res)
        .map { qry_id, map_m8, ref_id, anno_tsv ->
            tuple(qry_id, map_m8, anno_tsv)
        }

    res = assign_taxid_ko_to_contigs(joined_ch)

    emit:
    contig_anno = res.contig_anno
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

workflow ANNOTATE_CONTIG_ALL {
    p = createNullParamsChannel()

    // 単体実行時も本番同様に tuple(qry_id, f) の 2要素チャンネルを生成
    contig_map_ch = Channel.fromPath(params.annotate_contig_m8_path)
                            .map { f -> 
                                def qry_id = f.name.replaceAll(/_contig_map\.m8$/, '')
                                tuple(qry_id, f) 
                            }
    annotated_ch  = Channel.fromPath(params.annotate_contig_tsv_path)
                            .map { f -> tuple('merged_all_samples', f) }

    out_ch = ANNOTATE_CONTIG_SUB(p, contig_map_ch, annotated_ch)
    out_ch.contig_anno.view { i -> "CONTIG ANNO TSV: $i" }
}

workflow {
    ANNOTATE_CONTIG_ALL()
}