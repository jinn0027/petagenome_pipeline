#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.diamond_diamond_makerefdb_memory = params.memory
params.diamond_diamond_makerefdb_threads = params.threads

params.diamond_diamond_blastp_memory = params.memory
params.diamond_diamond_blastp_threads = params.threads

params.diamond_task = "megadiamond"
params.diamond_num_alignments = "1"
params.diamond_perc_identity = "95"
params.diamond_evalue = "1e-20"
params.diamond_outfmt = 6

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process diamond_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.diamond_diamond_makerefdb_memory}"
    def threads = "${params.diamond_diamond_makerefdb_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${ref_id}
        diamond \\
            makedb \\
            --threads ${threads} \\
            --in ${ref} \\
            -d ${ref_id}/ref
        """
}

process diamond_blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.diamond_diamond_blastp_memory}"
    def threads = "${params.diamond_diamond_blastp_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.tsv", arity: '1')
    script:
        processProfile(task)
        """
        mkdir -p ${qry_id}
        diamond \\
            blastp \\
            --threads ${threads} \\
            -q ${qry} \\
            -d ${ref_db}/ref \\
            -o ${qry_id}/out.tsv
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// DIAMOND DB (.dmnd) 作成処理の本体
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    // [ p_val, ref_id, ref_path ] にフラット化
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = diamond_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// DIAMOND 検索処理 (diamond blastp) の本体
workflow MAP_SUB {
    take:
    p
    ref_db
    qry

    main:
    // ref_db と qry を結合してフラット化
    in_ch = ref_db.combine(qry).map { ref_id, ref_path, qry_id, qry_path ->
        tuple(ref_id, ref_path, qry_id, qry_path)
    }

    // パラメータ p と結合してフラット化
    in_ch = p.combine(in_ch).map { p_val, ref_id, ref_path, qry_id, qry_path ->
        tuple(p_val, ref_id, ref_path, qry_id, qry_path)
    }

    out = diamond_blastp(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. DB 作成のみ実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.diamond_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// B. 作成済み DB を使って検索のみ実行 (-entry MAP_ONLY)
workflow MAP_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.diamond_db) // 既存の .dmnd パス
    qry    = createSeqsChannel(params.diamond_qry)

    out_ch = MAP_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// C. DB作成 〜 検索を一括実行 (デフォルト または -entry DIAMOND_ALL)
workflow DIAMOND_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.diamond_ref)
    qry = createSeqsChannel(params.diamond_qry)

    // DB を作成して検索処理へ流し込む
    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    DIAMOND_ALL()
}
