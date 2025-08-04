#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.prodigal_prodigal_memory = params.memory
params.prodigal_prodigal_threads = params.threads

params.prodigal_procedure = "meta"
//params.prodigal_procedure = "single"
params.prodigal_outfmt = "gbk"
//params.prodigal_outfmt = "gff"
//params.prodigal_outfmt = "sco"

process prodigal {
    tag "${read_id}"
    container = "${params.petagenomeDir}/modules/prodigal/prodigal.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.prodigal_prodigal_memory} GB"
    cpus "${params.prodigal_prodigal_threads}"
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("${read_id}/out.faa", arity: '1'), \
              path("${read_id}/out.fna", arity: '1'), \
              path("${read_id}/out.${params.prodigal_outfmt}", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        read_=${read}
        echo ${read} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            read_=\${read_%%.gz}
            unpigz -c ${read} > \${read_}
        fi
        mkdir -p ${read_id}
        prodigal \\
            -p ${params.prodigal_procedure} \\
            -i \${read_} \\
            -f ${params.prodigal_outfmt} \\
            -a ${read_id}/out.faa \\
            -d ${read_id}/out.fna \\
            -o ${read_id}/out.${params.prodigal_outfmt}
        """
}

workflow {
    read = createSeqsChannel(params.test_prodigal_read)
    out = prodigal(read)
    out.view { i -> "${i}" }
}
