#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metaphlan_metaphlan_memory = params.memory
params.metaphlan_metaphlan_threads = params.threads

params.metaphlan_input_type = "fastq"
//fastq,fasta,bowtie2out,sam
params.metaphlan_db = "/dev/shm/petagenome_pipeline/external/metaphlan_db"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process metaphlan {
    tag "${read_id}"
    def local_db = "/opt/db"
    container = "${params.petagenomeDir}/modules/metaphlan/metaphlan.sif"
    containerOptions = "${params.apptainerRunOptions} -B ${params.metaphlan_db}:${local_db}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.metaphlan_metaphlan_memory}"
    def threads = "${params.metaphlan_metaphlan_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("${read_id}/out.sam", arity: '1'), \
              path("${read_id}/out.prof", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${read_id}
        metaphlan \\
            --nproc ${threads} \\
            --bowtie2db ${local_db} \\
            --input_type ${getParam(p, 'metaphlan_input_type')} \\
            --bowtie2out ${read_id}/out.sam \\
            ${read} ${read_id}/out.prof
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// MetaPhlAn 処理の本体
workflow METAPHLAN_SUB {
    take:
    p
    read

    main:
    // [ p_val, seq_id, seq_path ] にフラット化
    in_ch = p.combine(read).map { p_val, seq_id, seq_path ->
        tuple(p_val, seq_id, seq_path)
    }

    out = metaphlan(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow METAPHLAN_ALL {
    p    = createNullParamsChannel()
    read = createSeqsChannel(params.metaphlan_read)

    out_ch = METAPHLAN_SUB(p, read)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    METAPHLAN_ALL()
}
