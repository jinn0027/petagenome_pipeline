#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.error_correction_get_length_memory = params.memory
params.error_correction_get_length_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"
include { spades_error_correction } from "${params.petagenomeDir}/nf/lv1/spades"
include { fastqc } from "${params.petagenomeDir}/nf/lv1/fastqc"

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.error_correction_get_length_memory}"
    def threads = "${params.error_correction_get_length_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(id), path(reads, arity: '1..*')
    output:
        tuple val(id), path("${id}/*.length.txt")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${id}
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' \${i} | \
            python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fastq \
            > ${id}/\${i}.length.txt
        done
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// エラー修正 ＋ QC ＋ 配列長取得の一連処理
workflow ERROR_CORRECTION_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    // 1. SPAdes によるエラー修正
    ec_out = spades_error_correction(in_ch)

    // 修正済みペアリードの抽出 [ p_val, pair_id, paired_reads ]
    ec_paired_ch = p.combine(
        ec_out.map { pair_id, paired_reads, unpaired_reads ->
            tuple(pair_id, paired_reads)
        }
    ).map { p_val, pair_id, paired_reads ->
        tuple(p_val, pair_id, paired_reads)
    }

    // 2. FastQC と 配列長取得 (get_length)
    fqc_out = fastqc(ec_paired_ch)
    len_out = get_length(ec_paired_ch)

    emit:
    ec  = ec_out
    fqc = fqc_out
    len = len_out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow ERROR_CORRECTION_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.error_correction_reads)

    out_ch = ERROR_CORRECTION_SUB(p, reads)

    // 各出力チャンネルの確認
    out_ch.ec.view  { i -> "EC: $i" }
    out_ch.fqc.view { i -> "FQC: $i" }
    out_ch.len.view { i -> "LEN: $i" }
}

// デフォルトエントリーポイント
workflow {
    ERROR_CORRECTION_ALL()
}
