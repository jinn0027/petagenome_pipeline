#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.diamond_diamond_makerefdb_memory = params.memory
params.diamond_diamond_makerefdb_threads = params.threads

params.diamond_diamond_blastp_memory = params.memory
params.diamond_diamond_blastp_threads = params.threads

params.diamond_task = "megadiamond"
params.diamond_num_alignments = "1"
params.diamond_perc_identity = "95"
params.diamond_evalue = "1e-20"
params.diamond_outfmt = 6

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process diamond_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.diamond_diamond_makerefdb_memory}"
    def threads = "${params.diamond_diamond_makerefdb_threads}"
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
        mkdir -p ${ref_id}
        diamond \\
            makedb \\
            --threads ${threads} \\
            --in ${ref} \\
            -d ${ref_id}/ref
        """
}

process diamond_blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.diamond_diamond_blastp_memory}"
    def threads = "${params.diamond_diamond_blastp_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.tsv", arity: '1')
    script:
        processProfile(task)
        """
        mkdir -p ${qry_id}
        diamond \\
            blastp \\
            --threads ${threads} \\
            -q ${qry} \\
            -d ${ref_db}/ref \\
            -o ${qry_id}/out.tsv
        """
}

workflow {
    p = createNullParamsChannel()
    ref = createSeqsChannel(params.test_diamond_ref)
    qry = createSeqsChannel(params.test_diamond_qry)

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = diamond_makerefdb(p.combine(ref))
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = diamond_blastp(p.combine(in))
    out.view { i -> "$i" }
}
