#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bwa_makedb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bwa/bwa.sif"
    publishDir "${params.output}/bwa/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')

    output:
        tuple val(ref_id), path("db")
    script:
        """
        mkdir -p db
        bwa \\
            index \\
            -p db/${ref_id} \\
            ${ref}
        """
}

process bwa_mem {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bwa/bwa.sif"
    publishDir "${params.output}/bwa/${ref_id}/${qry_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out/*.sam", arity: '1')
    script:
        """
        mkdir -p out
        bwa \\
            mem \\
            -t ${params.threads} \\
            ${db}/${ref_id} \\
            ${qry} \\
            > out/${ref_id}_@_${qry_id}.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_bwa_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_bwa_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    db = bwa_makedb(ref)
    //db.view { i -> "$i" }
    in = db.combine(qry)
    //in.view { i -> "$i" }
    out = bwa_mem(in)
    out.view { i -> "$i" }
}

