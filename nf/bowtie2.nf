#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bowtie2_makedb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/bowtie2/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')

    output:
        tuple val(ref_id), path("db")
    script:
        """
        mkdir -p db
        bowtie2-build \\
            --threads ${params.threads} \\
            --seed ${params.random_seed} \\
            ${ref} \\
            db/${ref_id}
        """
}

process bowtie2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/bowtie2/${ref_id}/${qry_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out/out.sam")
    script:
        """
        mkdir -p out
        bowtie2 \\
            -p ${params.threads} \\
            --seed ${params.random_seed} \\
            -f \\
            -x ${db}/${ref_id} \\
            -U ${qry} \\
            > out/out.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_bowtie2_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_bowtie2_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    db = bowtie2_makedb(ref)
    //db.view { i -> "$i" }
    in = db.combine(qry)
    //in.view { i -> "$i" }
    out = bowtie2(in)
    //out.view { i -> "$i" }
}

