#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.prinseq_trim_right = 10
params.prinseq_trim_left = 10
params.prinseq_qual_right = 20
params.prinseq_qual_left = 20
params.prinseq_qual_window = 20
params.prinseq_min_len = 75
params.prinseq_derep = 1
params.prinseq_lc_method = "dust"
params.prinseq_lc_threshold = 7
params.prinseq_trim_ns_right = 1
params.prinseq_ns_max_n = 0
params.prinseq_out_good = "good"
params.prinseq_out_bad = "bad"

process prinseq {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/prinseq/prinseq.sif"
    publishDir "${params.output}/prinseq/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("good*"), path("bad*")
    script:
        """
        read0=${reads[0]}
        read1=${reads[1]}

        echo ${reads[0]} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            read0=\${read0%%.gz}
            read1=\${read1%%.gz}
            unpigz -c ${reads[0]} > \${read0}
            unpigz -c ${reads[1]} > \${read1}
        fi

        prinseq-lite.pl \\
            -trim_left ${params.prinseq_trim_left} \\
            -trim_right ${params.prinseq_trim_right} \\
            -trim_qual_left ${params.prinseq_qual_left} \\
            -trim_qual_right ${params.prinseq_qual_right} \\
            -trim_qual_window ${params.prinseq_qual_window} \\
            -min_len ${params.prinseq_min_len} \\
            -derep ${params.prinseq_derep} \\
            -lc_method ${params.prinseq_lc_method} \\
            -lc_threshold ${params.prinseq_lc_threshold} \\
            -trim_ns_right ${params.prinseq_trim_ns_right} \\
            -ns_max_n ${params.prinseq_ns_max_n} \\
            -out_good ${params.prinseq_out_good} \\
            -out_bad ${params.prinseq_out_bad} \\
            -fastq \${read0} \\
            -fastq2 \${read1}
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_prinseq_reads, checkIfExists: true)
    prinseq = prinseq(reads)
    //prinseq.view { i -> "${i}" }
}
