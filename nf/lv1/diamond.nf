#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. DIAMOND 固有の上限値定義
def DIAMOND_MAKEREFDB_MAX_MEMORY  = 32 // GB (バイナリインデックス構築用)
def DIAMOND_MAKEREFDB_MAX_THREADS = 8  // インデックス構築の並列化飽和点

def DIAMOND_SEARCH_MAX_MEMORY     = 64 // GB (ブロックキャッシュ・高速バッチ検索用)
def DIAMOND_SEARCH_MAX_THREADS    = 16 // メモリ帯域とスレッド並列効率の上限

// 3. 上限値による動的クリッピング
params.diamond_diamond_makerefdb_memory  = Math.min(params.memory as Integer, DIAMOND_MAKEREFDB_MAX_MEMORY)
params.diamond_diamond_makerefdb_threads = Math.min(params.threads as Integer, DIAMOND_MAKEREFDB_MAX_THREADS)

params.diamond_diamond_blastp_memory     = Math.min(params.memory as Integer, DIAMOND_SEARCH_MAX_MEMORY)
params.diamond_diamond_blastp_threads    = Math.min(params.threads as Integer, DIAMOND_SEARCH_MAX_THREADS)

params.diamond_diamond_blastx_memory     = Math.min(params.memory as Integer, DIAMOND_SEARCH_MAX_MEMORY)
params.diamond_diamond_blastx_threads    = Math.min(params.threads as Integer, DIAMOND_SEARCH_MAX_THREADS)

params.diamond_task = "megadiamond"
params.diamond_num_alignments = "1"
params.diamond_perc_identity = "95"
params.diamond_evalue = "1e-20"
params.diamond_outfmt = 6

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

// ------------------------------------------------------------------
// Processes
// ------------------------------------------------------------------

process diamond_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
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

// Protein vs Protein Search (blastp)
process diamond_blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output
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

// Nucleotide vs Protein Search (blastx - 6-frame translation)
process diamond_blastx {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.diamond_diamond_blastx_memory}"
    def threads = "${params.diamond_diamond_blastx_threads}"
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
            blastx \\
            --threads ${threads} \\
            -q ${qry} \\
            -d ${ref_db}/ref \\
            -o ${qry_id}/out.tsv
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// DIAMOND DB (.dmnd) 作成処理（常にタンパク質FASTAが対象）
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = diamond_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// タンパク質クエリ用 (diamond blastp)
workflow MAP_PROT_SUB {
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

    out = diamond_blastp(in_ch)

    emit:
    out = out
}

// 塩基クエリ用 (diamond blastx)
workflow MAP_NUCL_SUB {
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

    out = diamond_blastx(in_ch)

    emit:
    out = out
}

// 後方互換性維持のためのエイリアス (MAP_SUB -> MAP_PROT_SUB)
workflow MAP_SUB {
    take: p; ref_db; qry
    main: out = MAP_PROT_SUB(p, ref_db, qry)
    emit: out = out.out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// DB 作成のみ実行
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.diamond_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// タンパク質クエリでの検索のみ (既存 DB 利用)
workflow MAP_PROT_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.diamond_db)
    qry    = createSeqsChannel(params.diamond_qry)

    out_ch = MAP_PROT_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// 塩基クエリでの検索のみ (既存 DB 利用)
workflow MAP_NUCL_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.diamond_db)
    qry    = createSeqsChannel(params.diamond_qry)

    out_ch = MAP_NUCL_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// DB作成 〜 タンパク質検索 一括実行
workflow MAP_PROT_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.diamond_ref)
    qry = createSeqsChannel(params.diamond_qry)

    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_PROT_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// DB作成 〜 塩基検索 一括実行
workflow MAP_NUCL_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.diamond_ref)
    qry = createSeqsChannel(params.diamond_qry)

    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_NUCL_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント (タンパク質アライメントをデフォルトに設定)
workflow {
    MAP_PROT_ALL()
}