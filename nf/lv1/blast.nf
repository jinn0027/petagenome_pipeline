#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.blast_blast_makerefdb_memory = params.memory
params.blast_blast_makerefdb_threads = params.threads

params.blast_blastn_memory = params.memory
params.blast_blastn_threads = params.threads

params.blast_dbtype = "nucl"
params.blast_task = "megablast"
params.blast_num_alignments = "1"
params.blast_perc_identity = "95"
params.blast_evalue = "1e-20"
params.blast_strand = "both"
params.blast_outfmt = 6

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process blast_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.blast_blast_makerefdb_memory}"
    def threads = "${params.blast_blast_makerefdb_threads}"
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
        ref_=${ref}
        echo ${ref} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            ref_=\${ref_%%.gz}
            unpigz -c ${ref} > \${ref_}
        fi
        mkdir -p ${ref_id}
        ref_siz=\$( du -b \$(readlink -f \${ref_}) | awk '{print(\$1)}' )
        if [ \${ref_siz} -gt 0 ] ; then
            makeblastdb \\
                -in \${ref_} \\
                -out ${ref_id}/ref \\
                -dbtype ${getParam(p, 'blast_dbtype')} \\
                -parse_seqids
        fi
        """
}

process blastn {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.blast_blastn_memory}"
    def threads = "${params.blast_blastn_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        qry_=${qry}
        echo ${qry} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            qry_=\${qry_%%.gz}
            unpigz -c ${qry} > \${qry_}
        fi
        mkdir -p ${qry_id}
        touch ${qry_id}/out
        qry_siz=\$( du -b \$(readlink -f \${qry_}) | awk '{print(\$1)}' )
        db_siz=\$(ls ${ref_db} 2>/dev/null | wc -l)
        if [ \${qry_siz} -gt 0 ] && [ \${db_siz} -gt 0 ]; then
            blastn \\
                -task ${getParam(p, 'blast_task')} \\
                -num_threads ${threads} \\
                -query \${qry_} \\
                -db ${ref_db}/ref \\
                -perc_identity ${getParam(p, 'blast_perc_identity')} \\
                -evalue ${getParam(p, 'blast_evalue')} \\
                -outfmt ${getParam(p, 'blast_outfmt')} \\
                -num_alignments ${getParam(p, 'blast_num_alignments')} \\
                -strand  ${getParam(p, 'blast_strand')} \\
                -out ${qry_id}/out
        fi
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// DB作成処理の本体
workflow BUILD_DB_SUB {
    take:
    p
    ref

    main:
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = blast_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// BLAST検索処理の本体
workflow SEARCH_SUB {
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

    out = blastn(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン（-entry）用エントリーポイント
//    ※ take: を持たせず、params からチャンネルを生成する
// ==========================================

// A. DB作成のみ実行 (-entry BUILD_DB_ONLY)
workflow BUILD_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.test_blast_ref)

    BUILD_DB_SUB(p, ref)
}

// B. 作成済みDBで検索のみ実行 (-entry SEARCH_ONLY)
workflow SEARCH_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.test_blast_db) // 既存DBのパス指定
    qry    = createSeqsChannel(params.test_blast_qry)

    out_ch = SEARCH_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

// C. 一括実行 (デフォルト または -entry BLAST_ALL)
workflow BLAST_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.test_blast_ref)
    qry = createSeqsChannel(params.test_blast_qry)

    // DBを作成して、その出力を検索処理へ流し込む
    db_ch  = BUILD_DB_SUB(p, ref)
    out_ch = SEARCH_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    BLAST_ALL()
}
