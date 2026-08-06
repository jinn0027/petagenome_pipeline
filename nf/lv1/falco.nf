#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.falco_falco_memory = params.memory
params.falco_falco_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process falco {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/falco/falco.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.falco_falco_memory}"
    def threads = "${params.falco_falco_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        falco \\
            -t ${threads} \\
            -o ${pair_id} \\
            ${reads[0]} \\
            ${reads[1]}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// Falco (FastQC) 処理の本体
workflow FALCO_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = falco(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow FALCO_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.falco_reads)

    out_ch = FALCO_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    FALCO_ALL()
}
// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// Falco (FastQC) 処理の本体
workflow FALCO_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = falco(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow FALCO_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.falco_reads)

    out_ch = FALCO_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    FALCO_ALL()
}
