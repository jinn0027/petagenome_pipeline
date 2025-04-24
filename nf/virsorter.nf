#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.virsorter_db_type = "refseq"
//params.virsorter_db = "virome"
params.virsorter_aligner = "blast"
//params.virsorter_aligner = "diamond"

process virsorter {
    tag "${read_id}"
    def local_db = "/opt/VirSorter/virsorter-data"
    def local_mga = "/opt/VirSorter/mga_linux_ia64"
    container = "${params.petagenomeDir}/modules/virsorter/virsorter.sif"
    containerOptions "-B ${params.virsorter_db}:${local_db} -B ${params.virsorter_mga}:${local_mga} --writable"
    publishDir "${params.output}/virsorter/${read_id}", mode: 'copy'
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), path("${params.virsorter_aligner}/out/*.csv", arity: '1')
    script:
        """
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
            --wdir ${params.virsorter_aligner}/out \\
            --data-dir ${local_db} \\
            --ncpu ${params.cpus} \\
            --fna \${read_}"
        if [ "${params.virsorter_aligner}" = "diamond" ] ; then
            args+=" --diamond"
        fi
        wrapper_phage_contigs_sorter_iPlant.pl \${args}
        """
}

workflow {
    read = channel.fromPath(params.test_virsorter_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    //read.view { i -> "$i" }
    out = virsorter(read)
    out.view { i -> "$i" }
}

