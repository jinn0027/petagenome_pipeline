#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================
// 0. グローバルフォールバックの設定
// ==========================================
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// アライナーパラメータの初期化
params.align_ref_aligner        = params.align_ref_aligner ?: "mmseqs2"
params.align_ref_is_prebuilt_db = params.align_ref_is_prebuilt_db ?: false

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. インポートパスの動的決定
// ==========================================

// PZLAST のリポジトリパスとスクリプト存在確認
def pz_script  = params.pzrepoDir ? "${params.pzrepoDir}/nf/lv1/pzlast.nf" : null
def use_pzlast = (params.align_ref_aligner == 'pzlast') && pz_script && file(pz_script).exists()

// 使用するアライナーのパスを決定 (pzlast の条件が揃わなければ mmseqs2.nf を選択)
def aligner_path = use_pzlast ? pz_script : "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// 塩基配列（ヌクレオチド）用のサブワークフローをインポート
include { BUILD_REF_DB_NUCL_SUB; MAP_NUCL_SUB } from "${aligner_path}"


// ==========================================
// 2. マッピングサブワークフロー（処理の本体）
// ==========================================

workflow ALIGN_CONTIGS_TO_ORFS_SUB {
    take:
    p
    target_ref_or_db  // マッピング先リファレンス (塩基配列 FASTA等) または ビルド済み DB インデックス
    reads             // クエリ側の塩基配列/リード [ qry_id, contigs.fasta 等 ]

    main:
    // --- A. DB (インデックス) の準備 ---
    if (params.align_ref_is_prebuilt_db) {
        target_db = target_ref_or_db
    } else {
        target_db = BUILD_REF_DB_NUCL_SUB(p, target_ref_or_db).ref_db
    }

    // --- B. 相同性検索 (fmt6/m8 形式) ---
    search_out = MAP_NUCL_SUB(p, target_db, reads) // 出力: [ ref_id, qry_id, out.m8 ]

    emit:
    out = search_out.out
}


// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

// A. DB インデックスの作成のみを実行
workflow BUILD_REF_DB_ONLY {
    p          = createNullParamsChannel()
    target_ref = createSeqsChannel(params.align_ref_fasta)

    db_out = BUILD_REF_DB_NUCL_SUB(p, target_ref)

    db_out.ref_db.view { id, db_path ->
        "[BUILD_REF_DB_ONLY] Created Index/DB (${params.align_ref_aligner}): ${id} -> ${db_path}"
    }
}

// B. 未構築FASTAから全行程を実行 (-entry ALIGN_CONTIGS_TO_ORFS_ALL)
workflow ALIGN_CONTIGS_TO_ORFS_ALL {
    p          = createNullParamsChannel()
    target_ref = createSeqsChannel(params.align_ref_fasta)
    reads      = createSeqsChannel(params.align_ref_reads)

    params.align_ref_is_prebuilt_db = false

    out_ch = ALIGN_CONTIGS_TO_ORFS_SUB(p, target_ref, reads)
    out_ch.out.view { i -> "ALIGNMENT OUT (m8): $i" }
}

// C. 作成済み DB を使用してマッピングを実行 (-entry ALIGN_CONTIGS_TO_ORFS_WITH_DB)
workflow ALIGN_CONTIGS_TO_ORFS_WITH_DB {
    p         = createNullParamsChannel()
    target_db = createSeqsChannel(params.align_ref_prebuilt_db)
    reads     = createSeqsChannel(params.align_ref_reads)

    params.align_ref_is_prebuilt_db = true

    out_ch = ALIGN_CONTIGS_TO_ORFS_SUB(p, target_db, reads)
    out_ch.out.view { i -> "ALIGNMENT OUT (m8 / PREBUILT DB): $i" }
}

// デフォルトエントリーポイント
workflow {
    if (params.align_ref_is_prebuilt_db) {
        ALIGN_CONTIGS_TO_ORFS_WITH_DB()
    } else {
        ALIGN_CONTIGS_TO_ORFS_ALL()
    }
}