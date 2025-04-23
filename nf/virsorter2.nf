#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process virsorter2 {
    tag "${qry_id}"
    def local_db = "/opt/VirSorter2/db"
    container = "${params.petagenomeDir}/modules/virsorter2/virsorter2.sif"
    containerOptions "-B ${params.virsorter2_db} --writable"
    publishDir "${params.output}/virsorter2/${qry_id}", mode: 'copy'
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id) , path("out/*.tsv")
    script:
        """
        virsorter \\
            config \\
            --init-source \\
            --db-dir=${params.virsorter2_db}
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
    out = virsorter2(qry)
    //out.view { i -> "$i" }
}

