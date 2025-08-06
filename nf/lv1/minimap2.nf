#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.minimap2_minimap2_makerefdb_memory = params.memory
params.minimap2_minimap2_makerefdb_threads = params.threads

params.minimap2_minimap2_memory = params.memory
params.minimap2_minimap2_threads = params.threads

params.minimap2_minimap2_e2e_memory = params.memory
params.minimap2_minimap2_e2e_threads = params.threads

params.minimap2_ambiguous = "random"
params.minimap2_minid = 0.95
params.minimap2_pairlen = 1500

params.test_minimap2_e2e = false

process minimap2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.minimap2_minimap2_makerefdb_memory}"
    def threads = "${params.minimap2_minimap2_makerefdb_threads}"
    memory "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${ref_id}
        minimap2 \\
            -t ${threads} \\
            -a ${ref} \\
            -d ${ref_id}/ref.idx
        """
}

process minimap2 {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    def gb = "${params.minimap2_minimap2_memory}"
    def threads = "${params.minimap2_minimap2_threads}"
    memory "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${threads} \\
            -a ${ref_db}/ref.idx \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

process minimap2_e2e {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/minimap2/minimap2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy'
    def gb = "${params.minimap2_minimap2_e2e_memory}"
    def threads = "${params.minimap2_minimap2_e2e_threads}"
    memory "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(ref_id), path(ref, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        minimap2 \\
            -t ${threads} \\
            -a ${ref} \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

workflow {
    ref = createSeqsChannel(params.test_minimap2_ref)
    qry = createSeqsChannel(params.test_minimap2_qry)

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

