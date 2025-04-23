#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.diamond_task = "megadiamond"
params.diamond_num_alignments = "1"
params.diamond_perc_identity = "95"
params.diamond_evalue = "1e-20"
params.diamond_outfmt = 6

process diamond_makedb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/diamond/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')

    output:
        tuple val(ref_id), path("db")
    script:
        """
        mkdir -p db
        diamond \\
            makedb \\
            --threads ${params.threads} \\
            --in ${ref} \\
            -d db/${ref_id}
        """
}

process diamond_blastp {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/diamond/diamond.sif"
    publishDir "${params.output}/diamond/${ref_id}/${qry_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out/*.tsv", arity: '1')
    script:
        """
        mkdir -p out
        diamond \\
            blastp \\
            -q ${qry} \\
            -d ${db}/${ref_id} \\
            -o out/${ref_id}_@_${qry_id}.tsv
        """
}

workflow {
    ref = channel.fromPath(params.test_diamond_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }
    qry = channel.fromPath(params.test_diamond_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    db = diamond_makedb(ref)
    //db.view { i -> "$i" }
    in = db.combine(qry)
    //in.view { i -> "$i" }
    out = diamond_blastp(in)
    out.view { i -> "$i" }
}
