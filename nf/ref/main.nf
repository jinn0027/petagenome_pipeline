#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "$baseDir/fastq/*_R{1,2}_001.fastq.gz"
params.output = "output"
params.threads = 4
params.fastp_cut_mean_quality = 15
params.fastp_reads_minlength = 15
params.fastp_memory = 10
params.spades_memory = 20

process fastp {
  executor 'sge'
  tag "${pair_id}"
  clusterOptions "-S /bin/bash -l s_vmem=${params.fastp_memory}G,mem_req=${params.fastp_memory}G -pe def_slot ${params.threads} -cwd"
  publishDir "${params.output}/01_fastp", mode: 'copy'
  input:
    tuple val(pair_id), path(reads)
  output:
    tuple val(pair_id), path("${pair_id}_R*_001.fastp.fastq.gz")
  script:
    """
    /usr/local/package/fastp/0.23.4/bin/fastp -w ${params.threads} -y \\
    -i ${reads[0]} -o ${pair_id}_R1_001.fastp.fastq.gz \\
    -I ${reads[1]} -O ${pair_id}_R2_001.fastp.fastq.gz \\
    --cut_front --cut_tail \\
    --cut_mean_quality ${params.fastp_cut_mean_quality} \\
    --length_required ${params.fastp_reads_minlength}
    """
}

process spades {
  executor 'sge'
  tag "${pair_id}"
  clusterOptions "-S /bin/bash -l s_vmem=${params.spades_memory}G,mem_req=${params.spades_memory}G -pe def_slot ${params.threads} -cwd"
  publishDir "${params.output}/03_spades_assembly", mode: 'copy'
  input:
    tuple val(pair_id), path(reads)
  output:
    tuple val(pair_id), path("${pair_id}/scaffolds.fasta")
  script:
    """
    /usr/local/package/spades/3.15.5/bin/spades.py \\
        --meta \\
        --threads ${params.threads} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${pair_id}
    """
}

workflow {
    raw_short_reads = channel.fromFilePairs(params.reads, checkIfExists: true)
    fastp = fastp(raw_short_reads)
    spades = spades(fastp)
}
