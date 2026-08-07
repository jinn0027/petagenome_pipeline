#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.blast_blast_makerefdb_memory = params.memory
params.blast_blast_makerefdb_threads = params.threads

params.blast_blastn_memory = params.memory
params.blast_blastn_threads = params.threads

params.blast_dbtype = "nucl"
params.blast_program = "blastn"
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
    containerOptions = "${params.apptainerRunOptions}"
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
    containerOptions = "${params.apptainerRunOptions}"
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
                -strand ${getParam(p, 'blast_strand')} \\
                -out ${qry_id}/out
        fi
        """
}

process blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    containerOptions = "${params.apptainerRunOptions}"
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
            blastp \\
                -num_threads ${threads} \\
                -query \${qry_} \\
                -db ${ref_db}/ref \\
                -perc_identity ${getParam(p, 'blast_perc_identity')} \\
                -evalue ${getParam(p, 'blast_evalue')} \\
                -outfmt ${getParam(p, 'blast_outfmt')} \\
                -num_alignments ${getParam(p, 'blast_num_alignments')} \\
                -out ${qry_id}/out
        fi
        """
}

// ==========================================
// 1. サブワークフロー（基幹ロジック）
// ==========================================

// DB作成処理の本体
workflow BUILD_REF_DB_SUB {
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
workflow MAP_SUB {
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

    // p 内の blast_program の値に応じてプロセス分岐
    prog = p.map { map -> getParam(map, 'blast_program') }.first()

    if (prog == "blastp") {
        out = blastp(in_ch)
    } else {
        out = blastn(in_ch)
    }

    emit:
    out = out
}

// ==========================================
// 2. 塩基配列用 / タンパク質用 サブワークフロー（型明示型）
// ==========================================

// --- 塩基配列用 (Nucleotide: nucl / blastn) ---
workflow BUILD_REF_DB_NUCL_SUB {
    take: p; ref
    main:
    p_mod = p.map { map -> (map ?: [:]) + [blast_dbtype: "nucl"] }
    out = BUILD_REF_DB_SUB(p_mod, ref)
    emit: ref_db = out.ref_db
}

workflow MAP_NUCL_SUB {
    take: p; ref_db; qry
    main:
    p_mod = p.map { map -> (map ?: [:]) + [blast_dbtype: "nucl", blast_program: "blastn"] }
    out = MAP_SUB(p_mod, ref_db, qry)
    emit: out = out.out
}

// --- タンパク質用 (Protein: prot / blastp) ---
workflow BUILD_REF_DB_PROT_SUB {
    take: p; ref
    main:
    p_mod = p.map { map -> (map ?: [:]) + [blast_dbtype: "prot"] }
    out = BUILD_REF_DB_SUB(p_mod, ref)
    emit: ref_db = out.ref_db
}

workflow MAP_PROT_SUB {
    take: p; ref_db; qry
    main:
    p_mod = p.map { map -> (map ?: [:]) + [blast_dbtype: "prot", blast_program: "blastp"] }
    out = MAP_SUB(p_mod, ref_db, qry)
    emit: out = out.out
}

// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// --- 汎用 (params に依存) ---
workflow BUILD_REF_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.blast_ref)

    BUILD_REF_DB_SUB(p, ref)
}

workflow MAP_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.blast_db)
    qry    = createSeqsChannel(params.blast_qry)

    out_ch = MAP_SUB(p, ref_db, qry)
    out_ch.out.view { i -> "$i" }
}

workflow BLAST_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.blast_ref)
    qry = createSeqsChannel(params.blast_qry)

    db_ch  = BUILD_REF_DB_SUB(p, ref)
    out_ch = MAP_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// --- 塩基配列用エントリーポイント ---
workflow BLAST_NUCL_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.blast_ref)
    qry = createSeqsChannel(params.blast_qry)

    db_ch  = BUILD_REF_DB_NUCL_SUB(p, ref)
    out_ch = MAP_NUCL_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// --- タンパク質用エントリーポイント ---
workflow BLAST_PROT_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.blast_ref)
    qry = createSeqsChannel(params.blast_qry)

    db_ch  = BUILD_REF_DB_PROT_SUB(p, ref)
    out_ch = MAP_PROT_SUB(p, db_ch.ref_db, qry)

    out_ch.out.view { i -> "$i" }
}

// デフォルトエントリーポイント
workflow {
    BLAST_ALL()
}