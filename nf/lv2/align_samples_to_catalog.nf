#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// アライナーパラメータの初期化
params.align_samples_to_catalog_aligner        = "mmseqs2"
params.align_samples_to_catalog_is_prebuilt_db = false

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. インポートパスの動的決定
// ==========================================

// PZLAST のリポジトリパスとスクリプト存在確認
def pz_script = params.containsKey('pzrepoDir') && params.pzrepoDir ? "${params.pzrepoDir}/nf/lv1/pzlast.nf" : null
def use_pzlast = (params.align_samples_to_catalog_aligner == 'pzlast') && pz_script && file(pz_script).exists()

// 使用するアライナーのパスを決定 (pzlast の条件が揃わなければ mmseqs2.nf を選択)
def aligner_path = use_pzlast ? pz_script : "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// 塩基配列（ヌクレオチド）用のサブワークフローをインポート
include { BUILD_REF_DB_NUCL_SUB; MAP_NUCL_SUB } from "${aligner_path}"

// ==========================================
// 2. マッピングサブワークフロー（処理の本体）
// ==========================================

workflow ALIGN_SAMPLES_TO_CATALOG_SUB {
    take:
    p
    target_ref_or_db  // マッピング先リファレンス (塩基配列 FASTA等) または ビルド済み DB インデックス
    reads             // クエリ側の塩基配列/リード [ qry_id, samples.fasta 等 ]

    main:
    // --- A. DB (インデックス) の準備 ---
    def is_prebuilt_db_closure = { getParam(p, params, 'align_samples_to_catalog_is_prebuilt_db' ) }
    if (is_prebuilt_db_closure()) {
        target_db = target_ref_or_db
    } else {
        target_db = BUILD_REF_DB_NUCL_SUB(p, target_ref_or_db).ref_db
    }

    // --- B. 相同性検索 (fmt6/m8 形式) ---
    search_out = MAP_NUCL_SUB(p, target_db, reads) // 出力: [ ref_id, qry_id, out.m8 ]


    // --- C. ref_id を除外し、[ qry_id, out.m8 ] の形に整形して出力する
    clean_out = search_out.out.map { ref_id, qry_id, m8_file -> tuple(qry_id, m8_file) }

    emit:
    out = clean_out
}

// ==========================================
// 3. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

// A. DB インデックスの作成のみを実行
workflow BUILD_REF_DB_ONLY {
    p          = createNullParamsChannel()
    target_ref = createSeqsChannel(params.align_samples_to_catalog_fasta)

    db_out = BUILD_REF_DB_NUCL_SUB(p, target_ref)

    db_out.ref_db.view { id, db_path ->
        "[BUILD_REF_DB_ONLY] Created Index/DB (${params.align_samples_to_catalog_aligner}): ${id} -> ${db_path}"
    }
}

// B. 未構築FASTAから全行程を実行 (-entry ALIGN_SAMPLES_TO_CATALOG_ALL)
workflow ALIGN_SAMPLES_TO_CATALOG_ALL {
    p          = createNullParamsChannel()
    target_ref = createSeqsChannel(params.align_samples_to_catalog_fasta)
    reads      = createSeqsChannel(params.align_samples_to_catalog_reads)

    params.align_samples_to_catalog_is_prebuilt_db = false

    out_ch = ALIGN_SAMPLES_TO_CATALOG_SUB(p, target_ref, reads)
    out_ch.out.view { i -> "ALIGNMENT OUT (m8): $i" }
}

// C. 作成済み DB を使用してマッピングを実行 (-entry ALIGN_SAMPLES_TO_CATALOG_WITH_DB)
workflow ALIGN_SAMPLES_TO_CATALOG_WITH_DB {
    p         = createNullParamsChannel()
    target_db = createSeqsChannel(params.align_samples_to_catalog_prebuilt_db)
    reads     = createSeqsChannel(params.align_samples_to_catalog_reads)

    params.align_samples_to_catalog_is_prebuilt_db = true

    out_ch = ALIGN_SAMPLES_TO_CATALOG_SUB(p, target_db, reads)
    out_ch.out.view { i -> "ALIGNMENT OUT (m8 / PREBUILT DB): $i" }
}

// デフォルトエントリーポイント
workflow {
    if (params.align_samples_to_catalog_is_prebuilt_db) {
        ALIGN_SAMMPLES_TO_CATALOG_WITH_DB()
    } else {
        ALIGN_SAMPLES_TO_CATALOG_ALL()
    }
}