#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process virsorter2 {
    tag "${read_id}"
    def local_db = "/opt/VirSorter2/db"
    container = "${params.petagenomeDir}/modules/virsorter2/virsorter2.sif"
    containerOptions "-B ${params.virsorter2_db}:${local_db} -B /tmp:/home"
    //containerOptions "-B ${params.virsorter2_db}:${local_db} --writable-tmpfs"
    publishDir "${params.output}/${task.process}/${read_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("out/final-viral-boundary.tsv"), path("out/final-viral-score.tsv")
    script:
        """
        virsorter \\
            config \\
            --init-source \\
            --db-dir=${local_db}
        virsorter \\
            run \\
            -j ${params.threads} \\
            -w out \\
            -i ${read}
        """
}

workflow {
    read = channel.fromPath(params.test_virsorter2_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    read.view { i -> "$i" }
    out = virsorter2(read)
    out.view { i -> "$i" }
}

