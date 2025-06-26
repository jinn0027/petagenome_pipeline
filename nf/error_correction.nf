#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { spades_error_correction } from "${params.petagenomeDir}/nf/spades"
include { fastqc } from "${params.petagenomeDir}/nf/fastqc"

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(reads, arity: '1..*')
    output:
        tuple val(id), path("out/*.length.txt.gz")
    script:
        """
        mkdir -p out
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            zcat \$i | \
                awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' | \
                python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fastq | \
                pigz -c >> out/\$(basename \$i).length.txt.gz
        done
        """
}

workflow error_correction {
  take:
    reads
  main:
    ec = spades_error_correction(reads)
    fqc = fastqc( ec.map { id, reads, unpaired -> tuple( id, reads ) } )
    len = get_length(ec)
  emit:
    ec
    fqc
    len
}

workflow {
    reads = channel.fromFilePairs(params.test_error_correction_reads, checkIfExists: true)
    out = error_correction(reads)
    out.ec.view{ i -> "$i" }
    out.fqc.view{ i -> "$i" }
    out.len.view{ i -> "$i" }
}
