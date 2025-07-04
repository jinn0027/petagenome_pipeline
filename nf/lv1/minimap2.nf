#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.minimap2_ambiguous = "random"
params.minimap2_minid = 0.95
params.minimap2_pairlen = 1500

params.test_minimap2_e2e = false

process minimap2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref, arity: '1')

    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        mkdir -p ${ref_id}
        minimap2 \\
            -t ${params.threads} \\
            -a ${ref} \\
            -d ${ref_id}/ref.idx
        """
}

process minimap2 {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')

    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${params.threads} \\
            -a ${ref_db}/ref.idx \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

process minimap2_e2e {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${params.threads} \\
            -a ${ref} \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_minimap2_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_minimap2_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    if (params.test_minimap2_e2e) {
        in = ref.combine(qry)
        out = minimap2_e2e(in)
        out.view { i -> "$i" }
    } else {
        ref_db = minimap2_makerefdb(ref)
        //ref_db.view { i -> "$i" }
        in = ref_db.combine(qry)
        //in.view { i -> "$i" }
        out = minimap2(in)
        out.view { i -> "$i" }
    }
}

