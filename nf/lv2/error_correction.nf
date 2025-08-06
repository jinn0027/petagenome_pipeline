#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createPairsChannel } from "${params.petagenomeDir}/nf/common/utils"
include { spades_error_correction } from "${params.petagenomeDir}/nf/lv1/spades"
include { fastqc } from "${params.petagenomeDir}/nf/lv1/fastqc"

params.error_correction_get_length_memory = params.memory
params.error_correction_get_length_threads = params.threads

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.error_correction_get_length_memory} GB"
    def threads = "${params.error_correction_get_length_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(reads, arity: '1..*')
    output:
        tuple val(id), path("${id}/*.length.txt")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' \${i} | \
            python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fastq \
            > ${id}/\${i}.length.txt
        done
        """
}

workflow error_correction {
  take:
    reads
  main:
    ec = spades_error_correction(reads)
    fqc = fastqc( ec.map { id, reads, unpaired -> tuple( id, reads ) } )
    len = get_length( ec.map { id, reads, unpaired -> tuple( id, reads ) } )
  emit:
    ec
    fqc
    len
}

workflow {
    reads = createPairsChannel(params.test_error_correction_reads)
    out = error_correction(reads)
    out.ec.view{ i -> "$i" }
    out.fqc.view{ i -> "$i" }
    out.len.view{ i -> "$i" }
}
