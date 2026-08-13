#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def EXTRACT_MAX_MEMORY  = 8
def EXTRACT_MAX_THREADS = 4

// 3. 上限値による動的クリッピング
params.remove_host_extract_memory  = Math.min(params.memory as Integer, EXTRACT_MAX_MEMORY)
params.remove_host_extract_threads = Math.min(params.threads as Integer, EXTRACT_MAX_THREADS)

params.remove_host_aligner          = "bwa_mem2"
params.remove_host_is_prebuilt_db   = false

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// 1. インポートパスの動的決定
// ==========================================

// params.pzrepoDir が安全に設定されており、かつ pzbwa.nf が存在するかチェック
def pz_script = (params.containsKey('pzrepoDir') && params.pzrepoDir) ? "${params.pzrepoDir}/nf/lv1/pzbwa.nf" : null
def use_pzbwa = (params.remove_host_aligner == 'pzbwa') && pz_script && file(pz_script).exists()

// 使用するアライナーのパスを決定 (pzbwa を使う条件が揃っていなければ bwa_mem2.nf を選択)
def aligner_path = use_pzbwa ? pz_script : "${params.petagenomeDir}/nf/lv1/bwa_mem2.nf"

// 決定したパスからインデックス作成サブワークフローとマッピングサブワークフローをインポート
include { BUILD_REF_DB_SUB; MAP_SUB } from "${aligner_path}"


// ==========================================
// 2. プロセス定義 (未マッピングペア抽出)
// ==========================================

process EXTRACT_UNMAPPED_READS {
    tag "${qry_id}"

    container = "${params.petagenomeDir}/modules/samtools/samtools.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.remove_host_extract_memory}"
    def threads = "${params.remove_host_extract_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(ref_id), val(qry_id), path(bam_or_sam)

    output:
        tuple val(qry_id), path("${qry_id}_host_removed*.fastq.gz"), emit: reads

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt

    # 最初の一行のSAM FLAG (第2フィールド) から Paired-end (FLAG & 1) かを判定
    IS_PAIRED=\$(samtools view ${bam_or_sam} | head -n 1 | awk '{if (and(\$2, 1)) print "paired"; else print "single"}')

    if [ "\$IS_PAIRED" = "paired" ]; then
        # --- ペアエンド用処理 (-f 12 : 両方のリードが未マッピング) ---
        samtools fastq -f 12 -@ ${threads} \\
            -1 ${qry_id}_host_removed_R1.fastq.gz \\
            -2 ${qry_id}_host_removed_R2.fastq.gz \\
            -0 /dev/null -s /dev/null \\
            ${bam_or_sam}
    else
        # --- シングルエンド用処理 (-f 4 : 未マッピングリードのみ抽出) ---
        samtools fastq -f 4 -@ ${threads} \\
            ${bam_or_sam} | bgzip -@ ${threads} -c > ${qry_id}_host_removed_single.fastq.gz
    fi
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
    // --- A. DB (インデックス) の準備 ---
    if (params.remove_host_is_prebuilt_db) {
        host_db = host_ref_or_db
    } else {
        // インポート元が動的に切り替わっているため、そのまま呼び出しOK
        host_db = BUILD_REF_DB_SUB(p, host_ref_or_db).ref_db
    }

    // --- B. マッピングの実行 ---
    map_out = MAP_SUB(p, host_db, reads)

    // --- C. 未マッピング（非ホスト）リードの抽出 ---
    cleaned = EXTRACT_UNMAPPED_READS(map_out.out)

    emit:
    reads = cleaned.reads
}


// ==========================================
// 4. テスト・単体実行用エントリーポイント (-entry)
// ==========================================

// A. DB (BWA/PZBWA インデックス) の作成のみを実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p        = createNullParamsChannel()
    host_ref = createSeqsChannel(params.remove_host_ref_fasta)

    db_out = BUILD_REF_DB_SUB(p, host_ref)

    db_out.ref_db.view { id, db_path ->
        "[BUILD_REF_DB_ONLY] Created Index/DB (${params.remove_host_aligner}): ${id} -> ${db_path}"
    }
}

// B. 未構築FASTAから全行程を実行 (-entry REMOVE_HOST_ALL)
workflow REMOVE_HOST_ALL {
    p        = createNullParamsChannel()
    host_ref = createSeqsChannel(params.remove_host_ref_fasta)
    reads    = createPairsChannel(params.remove_host_reads)

    params.remove_host_is_prebuilt_db = false

    out_ch = REMOVE_HOST_SUB(p, host_ref, reads)
    out_ch.reads.view { i -> "HOST REMOVED READS: $i" }
}

// C. 作成済み DB を使用してマッピング〜除去を実行 (-entry REMOVE_HOST_WITH_DB)
workflow REMOVE_HOST_WITH_DB {
    p       = createNullParamsChannel()
    host_db = createSeqsChannel(params.remove_host_prebuilt_db)
    reads   = createPairsChannel(params.remove_host_reads)

    params.remove_host_is_prebuilt_db = true

    out_ch = REMOVE_HOST_SUB(p, host_db, reads)
    out_ch.reads.view { i -> "HOST REMOVED READS (PREBUILT DB): $i" }
}

// デフォルトエントリーポイント
workflow {
    if (params.remove_host_is_prebuilt_db) {
        REMOVE_HOST_WITH_DB()
    } else {
        REMOVE_HOST_ALL()
    }
}