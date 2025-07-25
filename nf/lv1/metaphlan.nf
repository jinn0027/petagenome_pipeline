#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.metaphlan_metaphlan_memory = params.memory
params.metaphlan_metaphlan_threads = params.threads

params.metaphlan_input_type = "fastq"
//fastq,fasta,bowtie2out,sam
params.metaphlan_db = "/dev/shm/petagenome_pipeline/external/metaphlan_db"

process metaphlan {
    tag "${read_id}"
    def local_db = "/opt/db"
    container = "${params.petagenomeDir}/modules/metaphlan/metaphlan.sif"
    containerOptions "${params.apptainerRunOptions} -B ${params.metaphlan_db}:${local_db} -B /tmp:/home"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.metaphlan_metaphlan_memory} GB"
    cpus "${params.metaphlan_metaphlan_threads}"

    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("${read_id}/out.sam", arity: '1'), \
              path("${read_id}/out.prof", arity: '1')
    script:
        """
        mkdir -p ${read_id}
        metaphlan \\
            --nproc ${params.metaphlan_metaphlan_threads} \\
            --bowtie2db ${local_db} \\
            --input_type ${params.metaphlan_input_type} \\
            --bowtie2out ${read_id}/out.sam \\
            ${read} ${read_id}/out.prof
        """
}

workflow {
    read = createSeqsChannel(params.test_metaphlan_read)
    //read.view { i -> "${i}" }
    out = metaphlan(read)
    out.view { i -> "${i}" }
}
