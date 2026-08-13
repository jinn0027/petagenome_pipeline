#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

include { createNullParamsChannel; getParam } from "${params.petagenomeDir}/nf/common/utils"
include { FASTP_SUB }                         from "${params.petagenomeDir}/nf/lv1/fastp"
include { ERROR_CORRECTION_SUB }              from "${params.petagenomeDir}/nf/lv2/error_correction"
include { ASSEMBLY_SUB }                      from "${params.petagenomeDir}/nf/lv2/assembly"
include { POOL_CONTIGS_SUB }                  from "${params.petagenomeDir}/nf/lv2/pool_contigs"
include { CIRCULAR_CONTIGS_SUB }              from "${params.petagenomeDir}/nf/lv2/circular_contigs"

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

workflow BACTERIOME_PIPELINE_SUB {
    take:
    p
    reads
    l_thre

    main:
    // 1. fastp (QC/Trimming)
    fp = FASTP_SUB(p, reads)

    // 2. Error Correction
    ec = ERROR_CORRECTION_SUB(p, fp.out)

    // エラー修正済みリードの抽出 [ pair_id, paired ]
    ec_paired = ec.ec.map { pair_id, paired, unpaired ->
        tuple(pair_id, paired)
    }

    // 3. Assembly
    as_out = ASSEMBLY_SUB(p, ec_paired, l_thre)

    // 4. アセンブリ後のコンティグ収集 (Pool Contigs 用の前処理)
    // 全サンプルの Contig パスを集約し、最初のサンプル ID を代表キーとしてタプル化
    flt_all = as_out.flt
        .collect()
        .map { list_of_tuples ->
            def first_key   = list_of_tuples[0][0]
            def all_contigs = list_of_tuples.collect { key, contigs -> contigs }
            tuple(first_key, all_contigs)
        }

    // 5. Pool Contigs (複数サンプルの統合・重複除去)
    pc = POOL_CONTIGS_SUB(p, flt_all, l_thre)

    // 6. Circular Contigs (環状コンティグ検出)
    pc_flt_seqs = pc.flt.map { id, fa, name -> tuple(id, fa) }
    cc = CIRCULAR_CONTIGS_SUB(p, pc_flt_seqs)

    emit:
    reads  = reads
    fp     = fp.out
    ec_ec  = ec.ec
    ec_fqc = ec.fqc
    ec_len = ec.len
    as_flt = as_out.flt
    as_len = as_out.len
    as_sts = as_out.sts
    pc_db  = pc.blstdb
    cc_ded = cc.dedupl
}

// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow BACTERIOME_PIPELINE_ALL {
    p      = createNullParamsChannel()
    l_thre = params.bacteriome_pipeline_lthre

    // 入力ファイルのパースと Mix 処理
    def reads_list = params.bacteriome_pipeline_reads.split(';')

    def individual_channels = reads_list.collect { reads_path ->
        Channel.fromFilePairs(reads_path, checkIfExists: true)
    }

    // 全チャンネルの結合
    def reads_mixed = individual_channels.first()
    if (individual_channels.size() > 1) {
        individual_channels.tail().each { ch ->
            reads_mixed = reads_mixed.mix(ch)
        }
    }

    // mapWithIndex を使用した安全な連番 ID 化（並列実行時の競合回避）
    reads = reads_mixed.mapWithIndex { tuple_item, idx ->
        def id   = tuple_item[0]
        def pair = tuple_item[1]
        def new_id = "${String.format('%02d', idx)}_${id}"
        tuple(new_id, pair)
    }

    // パイプライン本体の呼び出し
    out_ch = BACTERIOME_PIPELINE_SUB(p, reads, l_thre)

    // 主要な出力チャンネルの確認
    out_ch.reads.view  { i -> "INPUT READS: $i" }
    out_ch.as_flt.view { i -> "ASSEMBLY FLT: $i" }
    out_ch.cc_ded.view { i -> "CIRCULAR DEDUPL: $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}