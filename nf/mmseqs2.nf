#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.mmseqs2_search_type = "3"

process mmseqs2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/mmseqs2", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("ref_db/${ref_id}")
    script:
        """
        mkdir -p ref_db/${ref_id}
        mmseqs createdb \\
            ${ref} \\
            ref_db/${ref_id}/${ref_id}
        """
}

process mmseqs2_makeqrydb {
    tag "${qry_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/mmseqs2", mode: 'copy'
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id), path("qry_db/${qry_id}")
    script:
        """
        mkdir -p qry_db/${qry_id}
        mmseqs createdb \\
            ${qry} \\
            qry_db/${qry_id}/${qry_id}
        """
}

process mmseqs2 {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/mmseqs2/${ref_id}/${qry_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref_db), val(qry_id), path(qry_db)
    output:
        tuple val(ref_id), path("${ref_id}_@_${qry_id}.*")
    script:
        """
        mkdir -p out tmp
        mmseqs search \\
            --threads ${params.threads} \\
            --search-type ${params.mmseqs2_search_type} \\
            ${qry_db}/${qry_id} \\
            ${ref_db}/${ref_id} \\
            ${ref_id}_@_${qry_id} \\
            tmp
        """
}

workflow {
    ref = channel.fromPath(params.test_mmseqs2_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_mmseqs2_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = mmseqs2_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    qry_db = mmseqs2_makeqrydb(qry)
    //qry_db.view { i -> "$i" }

    in = ref_db.combine(qry_db)
    //in.view { i -> "$i" }

    out = mmseqs2(in)
    //out.view { i -> "$i" }
}

