#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"
// bwa_mem2 は常にインポート
include { BUILD_DB_SUB as BWA_BUILD; MAP_SUB as BWA_MAP } from "${params.petagenomeDir}/nf/lv1/bwa_mem2.nf"
// pzbwa は params.pzrepoDir が指定されており、ファイルが存在する場合のみ読込 (optional)
include { BUILD_DB_SUB as PZ_BUILD;  MAP_SUB as PZ_MAP  } from "${params.pzrepoDir ?: '.'}/nf/lv1/pzbwa.nf" optional true

// ==========================================
// 2. プロセス定義 (未マッピングペア抽出)
// ==========================================
process EXTRACT_UNMAPPED_READS {
    tag "${qry_id}"

    input:
    tuple val(p_val), val(ref_id), val(qry_id), path(bam_or_sam)

    output:
    tuple val(qry_id), path("${qry_id}_host_removed_R*.fastq.gz"), emit: reads

    script:
    """
    # -f 12 : ペアの両リードがホストゲノムに未マッピングのもののみ抽出
    samtools fastq -f 12 \
        -1 ${qry_id}_host_removed_R1.fastq.gz \
        -2 ${qry_id}_host_removed_R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        ${bam_or_sam}
    """
}


// ==========================================
// 3. ホスト除去サブワークフロー（処理の本体）
// ==========================================
workflow REMOVE_HOST_SUB {
    take:
    p
    host_ref_or_db  // FASTA または ビルド済み DB インデックス
    reads           // クリーン済みリードペア [ pair_id, [R1, R2] ]

    main:
    // ツール指定と pzbwa 利用可否の判定
    def tool = params.host_removal_tool ?: 'bwa_mem2'
    def has_pzbwa = params.pzrepoDir && file("${params.pzrepoDir}/nf/lv1/pzbwa.nf").exists()

    // pzbwa が指定されたが使えない場合の安全対策（警告を出して bwa_mem2 に倒す）
    if (tool == 'pzbwa' && !has_pzbwa) {
        log.warn "WARNING: 'pzbwa' requested, but 'params.pzrepoDir' is not set or 'pzbwa.nf' missing. Falling back to 'bwa_mem2'."
    }

    // --- A. DB (インデックス) の準備 ---
    if (params.host_is_prebuilt_db) {
        host_db = host_ref_or_db
    } else {
        if (tool == 'pzbwa' && has_pzbwa) {
            host_db = PZ_BUILD(p, host_ref_or_db).ref_db
        } else {
            host_db = BWA_BUILD(p, host_ref_or_db).ref_db
        }
    }

    // --- B. マッピングの実行 ---
    if (tool == 'pzbwa' && has_pzbwa) {
        map_out = PZ_MAP(p, host_db, reads)
    } else {
        map_out = BWA_MAP(p, host_db, reads)
    }

    // --- C. 未マッピング（非ホスト）リードの抽出 ---
    cleaned = EXTRACT_UNMAPPED_READS(map_out.out)

    emit:
    reads = cleaned.reads
}


// ==========================================
// 4. テスト・単体実行用エントリーポイント
// ==========================================
workflow REMOVE_HOST_ALL {
    p        = createNullParamsChannel()
    host_ref = createSeqsChannel(params.test_host_ref_fasta)
    reads    = createPairsChannel(params.test_host_removal_reads)

    params.host_is_prebuilt_db = false

    out_ch = REMOVE_HOST_SUB(p, host_ref, reads)
    out_ch.reads.view { i -> "HOST REMOVED READS: $i" }
}

workflow REMOVE_HOST_WITH_DB {
    p       = createNullParamsChannel()
    host_db = createSeqsChannel(params.test_host_prebuilt_db)
    reads   = createPairsChannel(params.test_host_removal_reads)

    params.host_is_prebuilt_db = true

    out_ch = REMOVE_HOST_SUB(p, host_db, reads)
    out_ch.reads.view { i -> "HOST REMOVED READS (PREBUILT DB): $i" }
}

workflow {
    REMOVE_HOST_ALL()
}
