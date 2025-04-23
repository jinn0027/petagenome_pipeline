#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.metaphlan2_input_type = "fastq"
//fastq,fasta,multifasta,multifastq,bowtie2out,sam
params.metaphlan2_db = "/dev/shm/petagenome_pipeline/external/metaphlan2_db"

process metaphlan2 {
    tag "${read_id}"
    def local_db = "/opt/MetaPhlAn2/db_v20"
    container = "${params.petagenomeDir}/modules/metaphlan2/metaphlan2.sif"
    containerOptions "-B ${db}:${local_db}"
    publishDir "${params.output}/metaphlan2/${read_id}", mode: 'copy'
    input:
        tuple path(db, arity: '1'), val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("out")
    script:
        """
        mkdir -p out
        metaphlan2.py \\
            --nproc ${params.threads} \\
            --bowtie2db ${local_db}/mpa_v20_m200 \\
            --input_type ${params.metaphlan2_input_type} \\
            --bowtie2out out/out.sam \\
            ${read} out/out.prof
        """
}

workflow {
    db = channel.fromPath(params.metaphlan2_db, checkIfExists: true)
    read = channel.fromPath(params.test_metaphlan2_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    in = db.combine(read)
    out = metaphlan2(in)
    //out.view { i -> "${i}" }
}
