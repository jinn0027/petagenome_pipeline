#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bwa_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bwa/bwa.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        mkdir -p ${ref_id}
        bwa \\
            index \\
            -p ${ref_id}/ref \\
            ${ref}
        """
}

process bwa_mem {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bwa/bwa.sif"
    publishDir "${params.output}/${task.process}/${ref_id}/${qry_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out.sam", arity: '1')
    script:
        """
        mkdir -p out
        bwa \\
            mem \\
            -t ${params.threads} \\
            ${ref_db}/ref \\
            ${qry} \\
            > out.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_bwa_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_bwa_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
 
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = bwa_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = bwa_mem(in)
    out.view { i -> "$i" }
}

