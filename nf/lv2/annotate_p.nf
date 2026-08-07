#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.annotate_p_aligner = "mmseqs2"
params.annotate_p_is_prebuilt_db = false
params.annotate_p_annotate_orfs_memory = params.memory ?: 16
params.annotate_p_annotate_orfs_threads = 2

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// アライナー（mmseqs2 または pzlast）のパスを動的に決定してインポート
def pzlast_script = (params.containsKey('pzrepoDir') && params.pzrepoDir) ? "${params.pzrepoDir}/nf/lv1/pzlast.nf" : null
def use_pzlast = (params.annotate_p_aligner == 'pzlast') && pzlast_script && file(pzlast_script).exists()

def aligner_path = use_pzlast ? pzlast_script : "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// mmseqs2.nf / pzlast.nf で定義されているタンパク質用サブワークフローをインポート
include { BUILD_REF_DB_PROT_SUB; MAP_PROT_SUB } from "${aligner_path}"

// Python によるアノテーション結合プロセス
process ANNOTATE_ORFS {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output

    def gb = "${params.annotate_p_annotate_orfs_memory}"
    def threads = "${params.annotate_p_annotate_orfs_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
    tuple val(ref_id), val(qry_id), path(fmt6_result)
    path taxid_map  // ワークフロー側で channel.value() にして渡す
    path ko_map     // ワークフロー側で channel.value() にして渡す

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
    taxid_map     // uniprot_to_taxid.tsv (File path or Channel)
    ko_map        // uniprot_to_ko.tsv (File path or Channel)

    main:
    // A. DB の準備（タンパク質用サブワークフローを呼び出し）
    if (params.annotate_p_is_prebuilt_db) {
        db = ref_or_db
    } else {
        db = BUILD_REF_DB_PROT_SUB(p, ref_or_db).ref_db
    }

    // B. 相同性検索 (タンパク質用サブワークフローを呼び出し / 出力: [ ref_id, qry_id, out.m8 ])
    search_out = MAP_PROT_SUB(p, db, orfs)

    // C. Map ファイルを value チャネル化して全サンプルに再利用可能にする
    ch_taxid = (taxid_map instanceof Channel) ? taxid_map : Channel.value(taxid_map)
    ch_ko    = (ko_map instanceof Channel)    ? ko_map    : Channel.value(ko_map)

    // D. TaxID / KO の紐づけ
    annotated_out = ANNOTATE_ORFS(search_out, ch_taxid, ch_ko)

    emit:
    annotated = annotated_out.annotated
}

// ==========================================
// テスト・単体実行用エントリーポイント
// ==========================================

// A. DB (MMseqs2/PZLAST インデックス) の作成のみを実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p              = createNullParamsChannel()
    annotate_p_ref = createSeqsChannel(params.annotate_p_ref_fasta ?: params.ref_fasta)

    db_out = BUILD_REF_DB_PROT_SUB(p, annotate_p_ref)

    db_out.ref_db.view { id, db_path ->
        "[BUILD_REF_DB_ONLY] Created Index/DB (${params.annotate_p_aligner}): ${id} -> ${db_path}"
    }
}

// B. FASTA からリファレンス DB を構築して検索・アノテーションを行うワークフロー
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

// C. 事前構築済み DB を使用して検索・アノテーションを行うワークフロー
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