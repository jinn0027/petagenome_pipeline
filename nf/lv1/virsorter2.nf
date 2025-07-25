#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.virsorter2_virsorter2_memory = params.memory
params.virsorter2_virsorter2_threads = params.threads

process virsorter2 {
    tag "${read_id}"
    def local_db = "/opt/VirSorter2/db"
    container = "${params.petagenomeDir}/modules/virsorter2/virsorter2.sif"
    containerOptions "${params.apptainerRunOptions} -B ${params.virsorter2_db}:${local_db} -B /tmp:/home"
    //containerOptions "${params.apptainerRunOptions} -B ${params.virsorter2_db}:${local_db} --writable-tmpfs"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.virsorter2_virsorter2_memory} GB"
    cpus "${params.virsorter2_virsorter2_threads}"

    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("${read_id}/final-viral-boundary.tsv"), path("${read_id}/final-viral-score.tsv")
    script:
        """
        virsorter \\
            config \\
            --init-source \\
            --db-dir=${local_db}
        virsorter \\
            run \\
            -j ${params.virsorter2_virsorter2_threads} \\
            -w ${read_id} \\
            -i ${read}
        """
}

workflow {
    read = createSeqsChannel(params.test_virsorter2_read)
    read.view { i -> "$i" }
    out = virsorter2(read)
    out.view { i -> "$i" }
}

