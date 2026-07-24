#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.prodigal_prodigal_memory = params.memory
params.prodigal_prodigal_threads = params.threads

params.prodigal_procedure = "meta"
//params.prodigal_procedure = "single"
params.prodigal_outfmt = "gbk"
//params.prodigal_outfmt = "gff"
//params.prodigal_outfmt = "sco"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process prodigal {
    tag "${read_id}"
    container = "${params.petagenomeDir}/modules/prodigal/prodigal.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.prodigal_prodigal_memory}"
    def threads = "${params.prodigal_prodigal_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("${read_id}/out.faa", arity: '1'), \
              path("${read_id}/out.fna", arity: '1'), \
              path("${read_id}/out.${params.prodigal_outfmt}", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        read_=${read}
        echo ${read} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            read_=\${read_%%.gz}
            unpigz -c ${read} > \${read_}
        fi
        mkdir -p ${read_id}
        prodigal \\
            -p ${getParam(p, 'prodigal_procedure')} \\
            -i \${read_} \\
            -f ${getParam(p, 'prodigal_outfmt')} \\
            -a ${read_id}/out.faa \\
            -d ${read_id}/out.fna \\
            -o ${read_id}/out.${getParam(p, 'prodigal_outfmt')}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// Prodigal (遺伝子予測) 処理の本体
workflow PRODIGAL_SUB {
    take:
    p
    read

    main:
    // [ p_val, seq_id, seq_path ] にフラット化
    in_ch = p.combine(read).map { p_val, seq_id, seq_path ->
        tuple(p_val, seq_id, seq_path)
    }

    out = prodigal(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow PRODIGAL_ALL {
    p    = createNullParamsChannel()
    read = createSeqsChannel(params.test_prodigal_read)

    out_ch = PRODIGAL_SUB(p, read)
    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    PRODIGAL_ALL()
}
