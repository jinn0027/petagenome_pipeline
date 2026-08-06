#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// アライナー（mmseqs2 または pzlast）のパスを動的に決定してインポート
def pzlast_script = (params.containsKey('pzrepoDir') && params.pzrepoDir) ? "${params.pzrepoDir}/nf/lv1/pzlast.nf" : null
def use_pzlast = (params.annotate_p_aligner == 'pzlast') && pzlast_script && file(pzlast_script).exists()

def aligner_path = use_pzlast ? pzlast_script : "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"
include { BUILD_REF_DB_SUB; BUILD_QRY_DB_SUB; MAP_SUB } from "${aligner_path}"

// Python によるアノテーション結合プロセス
process ANNOTATE_ORFS {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output

    input:
    tuple val(ref_id), val(qry_id), path(fmt6_result)
    path taxid_map
    path ko_map

    output:
    tuple val(qry_id), path("${qry_id}_annotated.tsv"), emit: annotated

    script:
    """
    python3 - << 'EOF'

import sys

taxid_dict = {}
with open("${taxid_map}") as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 2:
            taxid_dict[parts[0]] = parts[1]

ko_dict = {}
with open("${ko_map}") as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 2:
            ko_dict[parts[0]] = parts[1]

output_file = "${qry_id}_annotated.tsv"
with open("${fmt6_result}") as fin, open(output_file, 'w') as fout:
    fout.write("qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\ttaxid\\tko\\n")
    
    for line in fin:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\\t')
        if len(cols) < 2:
            continue
        
        sseqid = cols[1]
        taxid = taxid_dict.get(sseqid, "N/A")
        ko = ko_dict.get(sseqid, "N/A")
        
        fout.write(f"{line.strip()}\\t{taxid}\\t{ko}\\n")
EOF
    """
}

// アノテーション・サブワークフローの本体
workflow ANNOTATE_TAXID_KO_SUB {
    take:
    p
    ref_or_db     // リファレンスFASTA または ビルド済みDB
    orfs          // Prodigal等の出力から抽出した [ qry_id, out.faa ]
    taxid_map     // uniprot_to_taxid.tsv
    ko_map        // uniprot_to_ko.tsv

    main:
    // A. DB の準備
    if (params.annotate_p_is_prebuilt_db) {
        db = ref_or_db
    } else {
        db = BUILD_REF_DB_SUB(p, ref_or_db).ref_db
    }

    // B. クエリDBの作成
    qry_db = BUILD_QRY_DB_SUB(p, orfs)

    // C. 相同性検索 (fmt6/m8 形式)
    search_out = MAP_SUB(p, db, qry_db.qry_db) // 出力: [ ref_id, qry_id, out.m8 ]

    // D. TaxID / KO の紐づけ
    annotated_out = ANNOTATE_ORFS(search_out, taxid_map, ko_map)

    emit:
    annotated = annotated_out.annotated
}

// ==========================================
// テスト・単体実行用エントリーポイント
// ==========================================

// FASTA からリファレンス DB を構築して検索・アノテーションを行うワークフロー
workflow ANNOTATE_ALL {
    p         = createNullParamsChannel()
    ref_fasta = createSeqsChannel(params.annotate_p_ref_fasta)
    orfs      = createSeqsChannel(params.annotate_p_orfs)
    taxid_map = file(params.annotate_p_taxid_map)
    ko_map    = file(params.annotate_p_ko_map)

    params.annotate_p_is_prebuilt_db = false

    out_ch = ANNOTATE_TAXID_KO_SUB(p, ref_fasta, orfs, taxid_map, ko_map)
    out_ch.annotated.view { i -> "ANNOTATED RESULT: $i" }
}

// 事前構築済み DB を使用して検索・アノテーションを行うワークフロー
workflow ANNOTATE_WITH_DB {
    p         = createNullParamsChannel()
    ref_db    = createSeqsChannel(params.annotate_p_prebuilt_db)
    orfs      = createSeqsChannel(params.annotate_p_orfs)
    taxid_map = file(params.annotate_p_taxid_map)
    ko_map    = file(params.annotate_p_ko_map)

    params.annotate_p_is_prebuilt_db = true

    out_ch = ANNOTATE_TAXID_KO_SUB(p, ref_db, orfs, taxid_map, ko_map)
    out_ch.annotated.view { i -> "ANNOTATED RESULT (PREBUILT DB): $i" }
}

// メイン・デフォルトエントリーポイント
workflow {
    if (params.annotate_p_is_prebuilt_db) {
        ANNOTATE_WITH_DB()
    } else {
        ANNOTATE_ALL()
    }
}
