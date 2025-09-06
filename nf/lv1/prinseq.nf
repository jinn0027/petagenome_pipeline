#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.prinseq_prinseq_memory = params.memory
params.prinseq_prinseq_threads = params.threads

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

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process prinseq {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/prinseq/prinseq.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.prinseq_prinseq_memory}"
    def threads = "${params.prinseq_prinseq_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${pair_id}

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
            -out_good ${pair_id}/good \\
            -out_bad ${pair_id}/bad \\
            -fastq \${read0} \\
            -fastq2 \${read1}
        """
}

workflow {
    p = createNullParamsChannel()
    reads = createPairsChannel(params.test_prinseq_reads)
    out = prinseq(reads)
    out.view { i -> "${i}" }
}
