#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metaphlan_metaphlan_memory = params.memory
params.metaphlan_metaphlan_threads = params.threads

params.metaphlan_input_type = "fastq"
//fastq,fasta,bowtie2out,sam
params.metaphlan_db = "/dev/shm/petagenome_pipeline/external/metaphlan_db"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process metaphlan {
    tag "${read_id}"
    def local_db = "/opt/db"
    container = "${params.petagenomeDir}/modules/metaphlan/metaphlan.sif"
    containerOptions "${params.apptainerRunOptions} -B ${params.metaphlan_db}:${local_db} -B /tmp:/home"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.metaphlan_metaphlan_memory}"
    def threads = "${params.metaphlan_metaphlan_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("${read_id}/out.sam", arity: '1'), \
              path("${read_id}/out.prof", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${read_id}
        metaphlan \\
            --nproc ${threads} \\
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
