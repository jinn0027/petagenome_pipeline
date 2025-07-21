#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.bwa_bwa_makerefdb_memory = params.memory
params.bwa_bwa_makerefdb_threads = params.threads

params.bwa_bwa_mem_memory = params.memory
params.bwa_bwa_mem_threads = params.threads

process bwa_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bwa/bwa.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.bwa_bwa_makerefdb_memory} GB"
    cpus "${params.bwa_bwa_makerefdb_threads}"

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
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory "${params.bwa_bwa_mem_memory} GB"
    cpus "${params.bwa_bwa_mem_threads}"

    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        mkdir -p ${qry_id}
        bwa \\
            mem \\
            -t ${params.bwa_bwa_mem_threads} \\
            ${ref_db}/ref \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

workflow {
    ref = createSeqsChannel(params.test_bwa_ref)
    qry = createSeqsChannel(params.test_bwa_qry)
 
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = bwa_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = bwa_mem(in)
    out.view { i -> "$i" }
}

