#!/usr/init/env nextflow
nextflow.enable.dsl=2

// 1. デフォルトパラメータ
params.memory  = 16
params.threads = 4
params.output  = "./output"
params.publish_output = true
params.executor = "local"
params.apptainerRunOptions = ""

// ホスト除去用パラメータ
params.remove_host_ref_fasta_or_db = "/path/to/GRCh38/bwa_db"
params.remove_host_aligner   = "bwa_mem2"
params.remove_host_is_prebuilt_db = true

// 複数指定用のパラメータ（カンマ区切りで複数指定可能）
params.target_regions = "v1,v12,v34,v4"

// 2. モジュールのインポート
include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; \
          createSeqsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

include { FASTP_SUB       } from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB } from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { MEGAHIT_SUB     } from "${params.petagenomeDir}/nf/lv1/megahit.nf"
include { EXTRACT_16S_SUB } from "${params.petagenomeDir}/nf/lv1/barrnap.nf"

// ==========================================
// ワークフローの定義
// ==========================================
workflow WGS_TO_16S {
    take:
    p
    host_ref_or_db
    reads

    main:
    // 0. fastpによる品質管理・アダプタトリミング
    // 出力形式: [pair_id, [fastq_1, fastq_2]]
    filtered_reads = FASTP_SUB(p, reads)

    // 1. ホストゲノム除去 (fastp済みのリードを入力とする)
    host_removed_reads = REMOVE_HOST_SUB(p, host_ref_or_db, filtered_reads)

    // 2. MEGAHITでアセンブリ (出力: [id, contigs])
    asm_out = MEGAHIT_SUB(p, host_removed_reads.reads)

    // 3. ターゲット領域のチャンネルを作成 (例: ["v1", "v2", "v3", "v4"])
    regions_ch = Channel.from(params.target_regions.tokenize(','))
                        .map { it.trim() }

    // 4. アセンブリ結果と各領域を結合し、領域ごとに並行して16S抽出を実行
    asm_with_region = asm_out.out.combine(regions_ch)

    sim_16s_out = EXTRACT_16S_SUB(p, asm_with_region)

    emit:
    sim_16s_out = sim_16s_out
}

// ==========================================
// エントリーポイント
// ==========================================
workflow {
    p        = createNullParamsChannel()
    reads    = createPairsChannel(params.wgs_reads)
    host_ref = createSeqsChannel(params.remove_host_ref_fasta_or_db)

    out_ch = WGS_TO_16S(p, host_ref, reads)
    
    out_ch.sim_16s_out.view { id, region, fastqs -> 
        "Finished WGS-to-16S for $id [Region: $region] -> $fastqs" 
    }
}
