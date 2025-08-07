#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"

params.virsorter_virsorter_memory = params.memory
params.virsorter_virsorter_threads = params.threads

params.virsorter_db_type = "refseq"
//params.virsorter_db = "virome"
params.virsorter_aligner = "blast"
//params.virsorter_aligner = "diamond"

process virsorter {
    tag "${read_id}"
    def local_db = "/opt/VirSorter/virsorter-data"
    def local_mga = "/opt/VirSorter/mga_linux_ia64"
    container = "${params.petagenomeDir}/modules/virsorter/virsorter.sif"
    containerOptions "${params.apptainerRunOptions} -B ${params.virsorter_db}:${local_db} -B ${params.virsorter_mga}:${local_mga}"
    //containerOptions "--no-home -B ${params.virsorter_db}:${local_db} -B ${params.virsorter_mga}:${local_mga}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.virsorter_virsorter_memory}"
    def threads = "${params.virsorter_virsorter_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("${params.virsorter_aligner}/${read_id}/VIRSorter_global-phage-signal.csv", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        read_=${read}
        echo ${read} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            read_=\${read_%%.gz}
            unpigz -c ${read} > \${read_}
        fi
        db_type=${params.virsorter_db_type}
        if [ "\${db_type}" = "refseq" ] ; then
            db_type="1"
        elif [ "\${db_type}" = "virome" ] ; then
            db_type="2"
        fi
        args="\\
            --db \${db_type} \\
            --wdir ${params.virsorter_aligner}/${read_id} \\
            --data-dir ${local_db} \\
            --ncpu ${threads} \\
            --fna \${read_}"
        if [ "${params.virsorter_aligner}" = "diamond" ] ; then
            args+=" --diamond"
        fi
        wrapper_phage_contigs_sorter_iPlant.pl \${args}
        """
}

workflow {
    read = createSeqsChannel(params.test_virsorter_read)
    //read.view { i -> "$i" }
    out = virsorter(read)
    out.view { i -> "$i" }
}

