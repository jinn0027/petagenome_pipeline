#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bowtie2_bowtie2_makerefdb_memory = params.memory
params.bowtie2_bowtie2_makerefdb_threads = params.threads

params.bowtie2_bowtie2_memory = params.memory
params.bowtie2_bowtie2_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process bowtie2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.bowtie2_bowtie2_makerefdb_memory}"
    def threads = "${params.bowtie2_bowtie2_makerefdb_threads}"
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
        bowtie2-build \\
            --threads ${threads} \\
            --seed ${getParam(p, 'random_seed')} \\
            ${ref} \\
            ${ref_id}/ref
        """
}

process bowtie2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.bowtie2_bowtie2_memory}"
    def threads = "${params.bowtie2_bowtie2_threads}"
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
        bowtie2 \\
            -p ${threads} \\
            --seed ${getParam(p, 'random_seed')} \\
            -f \\
            -x ${ref_db}/ref \\
            -U ${qry} \\
            > ${qry_id}/out.sam
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// Bowtie2 リファレンス DB (Index) 作成処理の本体
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    // [ p_val, ref_id, ref_path ] にフラット化
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = bowtie2_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// Bowtie2 マッピング処理の本体
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

    out = bowtie2(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. DB 作成のみ実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.bowtie2_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// B. 作成済み DB を使ってマッピングのみ実行 (-entry MAP_ONLY)
workflow MAP_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.bowtie2_db) // 既存のインデックスパス
    qry    = createSeqsChannel(params.bowtie2_qry)

    out_ch = MAP_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// C. DB作成 〜 マッピングを一括実行 (デフォルト または -entry BOWTIE2_ALL)
workflow BOWTIE2_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.bowtie2_ref)
    qry = createSeqsChannel(params.bowtie2_qry)

    // インデックスを作成してマッピングへ流し込む
    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    BOWTIE2_ALL()
}
