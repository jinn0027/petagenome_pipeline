#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output = "output"
params.threads = 4
params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15
params.fastp_memory = 10
//params.spades_memory = 20

process fastp {
    tag "${pair_id}"
    publishDir "${params.output}/01_fastp", mode: 'copy'
    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path("${pair_id}_fastp_out*")
    script:
    """
    fastp -w ${params.threads} -y \\
    -i ${reads[0]} -o ${pair_id}_fastp_out1.fastq.gz \\
    -I ${reads[1]} -O ${pair_id}_fastp_out2.fastq.gz \\
    --cut_front --cut_tail \\
    --cut_mean_quality ${params.fastp_cut_mean_quality} \\
    --length_required ${params.fastp_reads_minlength}
    """
}

process test2 {
    publishDir "result2"
    
    input:
        tuple val(pair_id), path(infile)

    output:
        tuple val(pair_id), path("output2.txt")

    script:
    """
    echo "${infile}" > output2.txt
    
    """
}

workflow {
   raw_short_reads = channel.fromFilePairs(params.reads, checkIfExists: true)
   //raw_short_reads.view{ raw_short_pairs -> "$pairs"}
   ch_out = fastp(raw_short_reads)
   //ch_out.view{ out -> "$out"}
   //ch_out2 = test2(ch_out)
   //ch_out2.view{ out -> "$out"}
}
