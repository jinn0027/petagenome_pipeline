#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1コンティグ・1ORFあたりの保持上限数のデフォルト値
params.annotate_contig_top_n_contig = 3  // 1 Contig に対して保持する ORF の上限
params.annotate_contig_top_n_anno   = 3  // 1 ORF に対して保持する DB Hit (TaxID/KO) の上限

include { createNullParamsChannel; clusterOptions; processProfile } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. プロセス定義
// ==========================================

process ASSIGN_TAXID_KO_TO_CONTIGS_TOP_N {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.memory ?: 8}"
    def threads = "${params.cpus ?: 2}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
    tuple val(ref_id), val(qry_id), path(contig_map_m8), path(annotated_tsv)

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
        if (FNR == 1) next;
        orf_id   = \$1
        sseqid   = \$2
        pident   = \$3
        evalue   = \$11
        bitscore = \$12
        taxid    = \$13
        ko       = \$14

        if (anno_count[orf_id] < ANNO_TOP_N) {
            anno_count[orf_id]++
            idx = anno_count[orf_id]

            orf_sseqid[orf_id, idx]   = sseqid
            orf_pident[orf_id, idx]   = pident
            orf_evalue[orf_id, idx]   = evalue
            orf_bitscore[orf_id, idx] = bitscore
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
    contig_map_out // tuple(ref_id, qry_id, path_m8)
    annotated_res  // tuple(ref_id, path_tsv)

    main:
    // combine を使用して 1 (annotated_res) : N (contig_map_out) の全サンプルとの組み合わせを作る
    joined_ch = contig_map_out.combine(annotated_res)
        .map { map_tuple, anno_tuple -> 
            tuple(map_tuple[0], map_tuple[1], map_tuple[2], anno_tuple[1])
        }

    res = ASSIGN_TAXID_KO_TO_CONTIGS_TOP_N(joined_ch)

    emit:
    contig_anno = res.contig_anno
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

workflow ANNOTATE_CONTIG_ALL {
    p = createNullParamsChannel()

    contig_map_ch = Channel.fromPath(params.annotate_contig_m8_path)
                            .map { f -> 
                                def qry_id = f.name.replaceAll(/_contig_map\.m8\$/, '')
                                tuple('merged_all_samples', qry_id, f) 
                            }
    annotated_ch  = Channel.fromPath(params.annotate_contig_tsv_path)
                            .map { f -> tuple('merged_all_samples', f) }

    out_ch = ANNOTATE_CONTIG_SUB(p, contig_map_ch, annotated_ch)
    out_ch.contig_anno.view { i -> "CONTIG ANNO TSV: $i" }
}

workflow {
    ANNOTATE_CONTIG_ALL()
}