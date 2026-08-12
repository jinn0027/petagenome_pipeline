#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// 2. bbmap_makerefdb 固有の上限値定義
def BBMAP_MAKEREFDB_MAX_MEMORY  = 64 // GB (大型ゲノムインデックス作成を考慮)
def BBMAP_MAKEREFDB_MAX_THREADS = 16 // スレッド並列効率の頭打ちを考慮

// 3. 上限値による動的クリッピング
params.bbmap_bbmap_makerefdb_memory  = Math.min((params.bbmap_bbmap_makerefdb_memory  ?: params.memory) as Integer, BBMAP_MAKEREFDB_MAX_MEMORY)
params.bbmap_bbmap_makerefdb_threads = Math.min((params.bbmap_bbmap_makerefdb_threads ?: params.threads) as Integer, BBMAP_MAKEREFDB_MAX_THREADS)

params.bbmap_bbmap_memory = params.memory
params.bbmap_bbmap_threads = params.threads

params.bbmap_ambiguous = "random"
params.bbmap_minid = 0.95
params.bbmap_pairlen = 1500

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process bbmap_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    containerOptions = "= ${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.bbmap_bbmap_makerefdb_memory}"
    def threads = "${params.bbmap_bbmap_makerefdb_threads}"
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
        bbmap.sh \\
            -Xmx${gb}g \\
            threads=${threads} \\
            ref=${ref} \\
            path=${ref_id}
        """
}

process bbmap {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.bbmap_bbmap_memory}"
    def threads = "${params.bbmap_bbmap_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref_db, arity: '1'), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(ref_id), val(pair_id), path("${pair_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        bbmap.sh \\
            -Xmx${gb}g \\
            threads=${threads} \\
            ambiguous=${getParam(p, 'bbmap_ambiguous')} \\
            minid=${getParam(p, 'bbmap_minid')} \\
            pairlen=${getParam(p, 'bbmap_pairlen')} \\
            path=${ref_db} \\
            in=${reads[0]} \\
            in2=${reads[1]} \\
            out=${pair_id}/out.sam
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// BBMap リファレンス DB (Index) 作成処理の本体
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = bbmap_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// BBMap マッピング処理の本体
workflow MAP_SUB {
    take:
    p
    ref_db
    reads

    main:
    // ref_db と reads を結合してフラット化
    in_ch = ref_db.combine(reads).map { ref_id, ref_path, pair_id, reads_path ->
        tuple(ref_id, ref_path, pair_id, reads_path)
    }

    // パラメータ p と結合してフラット化
    in_ch = p.combine(in_ch).map { p_val, ref_id, ref_path, pair_id, reads_path ->
        tuple(p_val, ref_id, ref_path, pair_id, reads_path)
    }

    out = bbmap(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. DB 作成のみ実行 (-entry BUILD_REF_DB_ONLY)
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.bbmap_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// B. 作成済み DB を使ってマッピングのみ実行 (-entry MAP_ONLY)
workflow MAP_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.bbmap_db) // 既存のインデックスパス
    reads  = createPairsChannel(params.bbmap_reads)

    out_ch = MAP_SUB(p, ref_db, reads)
    out_ch.out.view { i -> "$i" }
}

// C. DB作成 〜 マッピングを一括実行 (デフォルト または -entry BBMAP_ALL)
workflow BBMAP_ALL {
    p     = createNullParamsChannel()
    ref   = createSeqsChannel(params.bbmap_ref)
    reads = createPairsChannel(params.bbmap_reads)

    // インデックスを作成してマッピングへ流し込む
    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_SUB(p, db_ch.ref_db, reads)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    BBMAP_ALL()
}
