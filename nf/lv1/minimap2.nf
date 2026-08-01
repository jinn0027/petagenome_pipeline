#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.minimap2_minimap2_makerefdb_memory = params.memory
params.minimap2_minimap2_makerefdb_threads = params.threads

params.minimap2_minimap2_memory = params.memory
params.minimap2_minimap2_threads = params.threads

params.minimap2_minimap2_e2e_memory = params.memory
params.minimap2_minimap2_e2e_threads = params.threads

params.minimap2_ambiguous = "random"
params.minimap2_minid = 0.95
params.minimap2_pairlen = 1500

params.minimap2_e2e = false

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process minimap2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.minimap2_minimap2_makerefdb_memory}"
    def threads = "${params.minimap2_minimap2_makerefdb_threads}"
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
        minimap2 \\
            -t ${threads} \\
            -a ${ref} \\
            -d ${ref_id}/ref.idx
        """
}

process minimap2 {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    def gb = "${params.minimap2_minimap2_memory}"
    def threads = "${params.minimap2_minimap2_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${threads} \\
            -a ${ref_db}/ref.idx \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

process minimap2_e2e {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    def gb = "${params.minimap2_minimap2_e2e_memory}"
    def threads = "${params.minimap2_minimap2_e2e_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${threads} \\
            -a ${ref} \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// Minimap2 リファレンス DB (.mmi インデックス) 作成処理の本体
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = minimap2_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// 作成済み DB (.mmi) を使った Minimap2 マッピング処理の本体
workflow MAP_SUB {
    take:
    p
    ref_db
    qry

    main:
    in_ch = ref_db.combine(qry).map { ref_id, ref_path, qry_id, qry_path ->
        tuple(ref_id, ref_path, qry_id, qry_path)
    }

    in_ch = p.combine(in_ch).map { p_val, ref_id, ref_path, qry_id, qry_path ->
        tuple(p_val, ref_id, ref_path, qry_id, qry_path)
    }

    out = minimap2(in_ch)

    emit:
    out = out
}

// DB作成なしで直接アライメントする End-to-End 処理の本体
workflow MAP_E2E_SUB {
    take:
    p
    ref
    qry

    main:
    in_ch = ref.combine(qry).map { ref_id, ref_path, qry_id, qry_path ->
        tuple(ref_id, ref_path, qry_id, qry_path)
    }

    in_ch = p.combine(in_ch).map { p_val, ref_id, ref_path, qry_id, qry_path ->
        tuple(p_val, ref_id, ref_path, qry_id, qry_path)
    }

    out = minimap2_e2e(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. DB (Index) 作成のみ実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.minimap2_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// B. 作成済み DB を使ってマッピングのみ実行 (-entry MAP_ONLY)
workflow MAP_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.minimap2_db) // 既存の .mmi パス
    qry    = createSeqsChannel(params.minimap2_qry)

    out_ch = MAP_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// C. 直接アライメントのみ実行 (-entry MAP_E2E_ONLY)
workflow MAP_E2E_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.minimap2_ref)
    qry = createSeqsChannel(params.minimap2_qry)

    out_ch = MAP_E2E_SUB(p, ref, qry)
    out_ch.out.view { i -> "$i" }
}

// D. フラグ分岐または一括実行 (デフォルト または -entry MINIMAP2_ALL)
workflow MINIMAP2_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.minimap2_ref)
    qry = createSeqsChannel(params.minimap2_qry)

    if (params.minimap2_e2e) {
        out_ch = MAP_E2E_SUB(p, ref, qry)
        out_ch.out.view { i -> "$i" }
    } else {
        db_ch  = BUILD_REF_DB_SUB(p, ref)
        out_ch = MAP_SUB(p, db_ch.ref_db, qry)
        out_ch.out.view { i -> "$i" }
    }
}

// デフォルトエントリーポイント
workflow {
    MINIMAP2_ALL()
}
