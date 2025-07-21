#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bowtie2_bowtie2_makerefdb_memory = params.memory
params.bowtie2_bowtie2_makerefdb_threads = params.threads

params.bowtie2_bowtie2_memory = params.memory
params.bowtie2_bowtie2_threads = params.threads

include { createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

process bowtie2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.bowtie2_bowtie2_makerefdb_memory} GB"
    cpus "${params.bowtie2_bowtie2_makerefdb_threads}"

    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        mkdir -p ${ref_id}
        bowtie2-build \\
            --threads ${params.bowtie2_bowtie2_makerefdb_threads} \\
            --seed ${params.random_seed} \\
            ${ref} \\
            ${ref_id}/ref
        """
}

process bowtie2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory "${params.bowtie2_bowtie2_memory} GB"
    cpus "${params.bowtie2_bowtie2_threads}"

    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        mkdir -p ${qry_id}
        bowtie2 \\
            -p ${params.bowtie2_bowtie2_threads} \\
            --seed ${params.random_seed} \\
            -f \\
            -x ${ref_db}/ref \\
            -U ${qry} \\
            > ${qry_id}/out.sam
        """
}

workflow {
    ref = createSeqsChannel(params.test_bowtie2_ref)
    qry = createSeqsChannel(params.test_bowtie2_qry)

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = bowtie2_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = bowtie2(in)
    out.view { i -> "$i" }
}

