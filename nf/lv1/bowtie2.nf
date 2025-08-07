#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.bowtie2_bowtie2_makerefdb_memory = params.memory
params.bowtie2_bowtie2_makerefdb_threads = params.threads

params.bowtie2_bowtie2_memory = params.memory
params.bowtie2_bowtie2_threads = params.threads

process bowtie2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.bowtie2_bowtie2_makerefdb_memory}"
    def threads = "${params.bowtie2_bowtie2_makerefdb_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${ref_id}
        bowtie2-build \\
            --threads ${threads} \\
            --seed ${params.random_seed} \\
            ${ref} \\
            ${ref_id}/ref
        """
}

process bowtie2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.bowtie2_bowtie2_memory}"
    def threads = "${params.bowtie2_bowtie2_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.sam", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        bowtie2 \\
            -p ${threads} \\
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

