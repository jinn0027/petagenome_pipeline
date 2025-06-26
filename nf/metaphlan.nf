#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metaphlan_input_type = "fastq"
//fastq,fasta,bowtie2out,sam
params.metaphlan_db = "/dev/shm/petagenome_pipeline/external/metaphlan2_db"

process metaphlan {
    tag "${read_id}"
    def local_db = "/opt/db"
    container = "${params.petagenomeDir}/modules/metaphlan/metaphlan.sif"
    containerOptions "-B ${params.metaphlan_db}:${local_db} -B /tmp:/home"
    publishDir "${params.output}/${task.process}/${read_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("out/*.sam", arity: '1'), \
              path("out/*.prof", arity: '1')
    script:
        """
        mkdir -p out
        metaphlan \\
            --nproc ${params.threads} \\
            --bowtie2db ${local_db} \\
            --input_type ${params.metaphlan_input_type} \\
            --bowtie2out out/${read_id}.sam \\
            ${read} out/${read_id}.prof
        """
}

workflow {
    read = channel.fromPath(params.test_metaphlan_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    //read.view { i -> "${i}" }
    out = metaphlan(read)
    out.view { i -> "${i}" }
}
