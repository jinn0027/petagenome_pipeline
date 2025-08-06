#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.bowtie_bowtie_makerefdb_memory = params.memory
params.bowtie_bowtie_makerefdb_threads = params.threads

params.bowtie_bowtie_memory = params.memory
params.bowtie_bowtie_threads = params.threads

process bowtie_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie/bowtie.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.bowtie_bowtie_makerefdb_memory} GB"
    threads = "${params.bowtie_bowtie_makerefdb_threads}"
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
        bowtie-build \\
            --threads ${params.bowtie_bowtie_makerefdb_threads} \\
            --seed ${params.random_seed} \\
            ${ref} \\
            ${ref_id}/ref
        """
}

process bowtie {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie/bowtie.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory "${params.bowtie_bowtie_memory} GB"
    threads = "${params.bowtie_bowtie_threads}"
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
        bowtie \\
            -p ${params.bowtie_bowtie_threads} \\
            --seed ${params.random_seed} \\
            -S \\
            -f \\
            -x ${ref_db}/ref \\
            ${qry} \\
            > ${qry_id}/out.sam
        """
}

workflow {
    ref = createSeqsChannel(params.test_bowtie_ref)
    qry = createSeqsChannel(params.test_bowtie_qry)

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = bowtie_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = bowtie(in)
    out.view { i -> "$i" }
}

