#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process virsorter2 {
    tag "${qry_id}"
    container = "${params.petagenomeDir}/modules/virsorter2/virsorter2.sif"
    publishDir "${params.output}/virsorter2/${qry_id}", mode: 'copy'
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id) , path("out/*.tsv")
    script:
        """
        virsorter \\
            run \\
            -j ${params.threads} \\
            -w out \\
            -i ${qry}
        """
}

workflow {
    qry = channel.fromPath(params.test_virsorter2_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
    qry.view { i -> "$i" }
    virsorter2 = virsorter2(qry)
    //virsorter2.view { i -> "$i" }
}

