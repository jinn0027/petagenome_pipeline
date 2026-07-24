#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createNullParamsChannel; createSeqsChannel; getParam } from "${params.petagenomeDir}/nf/common/utils"
include { FASTP_SUB }       from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB } from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { ASSEMBLY_SUB }    from "${params.petagenomeDir}/nf/lv2/assembly.nf"

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// バクテリオーム解析パイプライン
workflow BACTERIOME_PIPELINE_SUB {
    take:
    p
    host_ref  // ホストゲノム(FASTAまたはDB)のチャンネル
    reads     // 入力リードペア [ pair_id, [R1, R2] ]

    main:
    // 1. FASTP による QC・トリミング
    fp = FASTP_SUB(p, reads)

    // 2. ホスト除去 (FASTP通過後のクリーンリード fp.out を入力にする)
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)

    // 3. アセンブリ (長さの閾値をパラメータから取得。デフォルトは 1000 bp)
    l_thre = params.containsKey('assembly_l_thre') ? params.assembly_l_thre : 1000
    asm_res = ASSEMBLY_SUB(p, host_removed.reads, l_thre)

    emit:
    raw_reads   = reads
    fastp_reads = fp.out
    final_reads = host_removed.reads  // ホスト除去後の最終非ホストリード
    
    // アセンブリ結果の出力
    contigs     = asm_res.asm         // アセンブルされたコンティグ/スカッフォールド
    flt_contigs = asm_res.flt         // フィルタリング・リネーム後のコンティグ
    stats       = asm_res.sts         // 統計情報
    blastdb     = asm_res.blstdb      // 作成された BLAST DB
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow BACTERIOME_PIPELINE_ALL {
    p = createNullParamsChannel()

    // 1. 入力リードのチャンネル生成
    def reads_list = params.test_bacteriome_pipeline_reads.split(';')

    def individual_channels = reads_list.collect { reads_path ->
        channel.fromFilePairs(reads_path, checkIfExists: true)
    }

    def reads_mixed = individual_channels.first()
    individual_channels.tail().each { ch ->
        reads_mixed = reads_mixed.mix(ch)
    }

    def index = 0
    reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index++)}_${id}"
        tuple(new_id, pair)
    }

    // 2. ホスト参照配列（またはDB）のチャンネル生成
    host_ref = createSeqsChannel(params.host_ref_fasta)

    // 3. パイプライン本体の呼び出し
    out_ch = BACTERIOME_PIPELINE_SUB(p, host_ref, reads)

    // 結果の確認
    out_ch.fastp_reads.view { i -> "FASTP PASSED READS : $i" }
    out_ch.final_reads.view { i -> "HOST REMOVED READS : $i" }
    out_ch.contigs.view     { i -> "ASSEMBLY CONTIGS   : $i" }
    out_ch.stats.view       { i -> "ASSEMBLY STATS     : $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}
