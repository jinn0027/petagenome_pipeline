#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metaphlan_input_type = "fastq"
//fastq fasta bowtie2out sam

process metaphlan {
    tag "${read_id}"
    container = "${params.petagenomeDir}/modules/metaphlan/metaphlan.sif"
    publishDir "${params.output}/metaphlan/${read_id}", mode: 'copy'
    input:
        tuple path(db, arity: '1'), val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("out")
    script:
        """
        mkdir -p out
        metaphlan \\
            --nproc ${params.threads} \\
            --bowtie2db ${db} \\
            --input_type ${params.metaphlan_input_type} \\
            --bowtie2out out/out.sam \\
            ${read} out/out.prof
        """
}

workflow {
    db = channel.fromPath(params.test_metaphlan_db, checkIfExists: true)
    read = channel.fromPath(params.test_metaphlan_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    in = db.combine(read)
    out = metaphlan(in)
    //out.view { i -> "${i}" }
}
