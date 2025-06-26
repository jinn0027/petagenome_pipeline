#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bowtie2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif", enabled: params.publish_output
    publishDir "${params.output}/${task.process}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("ref_db/${ref_id}")
    script:
        """
        mkdir -p ref_db/${ref_id}
        bowtie2-build \\
            --threads ${params.threads} \\
            --seed ${params.random_seed} \\
            ${ref} \\
            ref_db/${ref_id}/${ref_id}
        """
}

process bowtie2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/bowtie2/bowtie2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}/${qry_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out/*.sam", arity: '1')
    script:
        """
        mkdir -p out
        bowtie2 \\
            -p ${params.threads} \\
            --seed ${params.random_seed} \\
            -f \\
            -x ${ref_db}/${ref_id} \\
            -U ${qry} \\
            > out/${ref_id}_${qry_id}.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_bowtie2_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_bowtie2_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = bowtie2_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = bowtie2(in)
    out.view { i -> "$i" }
}

