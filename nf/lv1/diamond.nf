#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.diamond_diamond_makerefdb_memory = params.memory
params.diamond_diamond_makerefdb_threads = params.threads

params.diamond_diamond_blastp_memory = params.memory
params.diamond_diamond_blastp_threads = params.threads

params.diamond_task = "megadiamond"
params.diamond_num_alignments = "1"
params.diamond_perc_identity = "95"
params.diamond_evalue = "1e-20"
params.diamond_outfmt = 6

process diamond_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.diamond_diamond_makerefdb_memory} GB"
    cpus "${params.diamond_diamond_makerefdb_threads}"
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${ref_id}
        diamond \\
            makedb \\
            --threads ${params.diamond_diamond_makerefdb_threads} \\
            --in ${ref} \\
            -d ${ref_id}/ref
        """
}

process diamond_blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory "${params.diamond_diamond_blastp_memory} GB"
    cpus "${params.diamond_diamond_blastp_threads}"

    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.tsv", arity: '1')
    script:
        processProfile(task)
        """
        mkdir -p ${qry_id}
        diamond \\
            blastp \\
            --threads ${params.diamond_diamond_blastp_threads} \\
            -q ${qry} \\
            -d ${ref_db}/ref \\
            -o ${qry_id}/out.tsv
        """
}

workflow {
    ref = createSeqsChannel(params.test_diamond_ref)
    qry = createSeqsChannel(params.test_diamond_qry)

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = diamond_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = diamond_blastp(in)
    out.view { i -> "$i" }
}
