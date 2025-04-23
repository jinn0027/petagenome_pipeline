#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.virsorter_db_type = "refseq"
//params.virsorter_db = "virome"
params.virsorter_aligner = "blast"
//params.virsorter_aligner = "diamond"

process virsorter {
    tag "${qry_id}"
    def local_db = "/opt/VirSorter/virsorter-data"
    def local_mga = "/opt/VirSorter/mga_linux_ia64"
    container = "${params.petagenomeDir}/modules/virsorter/virsorter.sif"
    containerOptions "-B ${params.virsorter_db}:${local_db} -B ${params.virsorter_mga}:${local_mga}"
    publishDir "${params.output}/virsorter/${qry_id}", mode: 'copy'
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id) , path("${params.virsorter_aligner}/out/*.csv", arity: '1')
    script:
        """
        qry_=${qry}
        echo ${qry} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            qry_=\${qry_%%.gz}
            unpigz -c ${qry} > \${qry_}
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
            --fna \${qry_}"
        if [ "${params.virsorter_aligner}" = "diamond" ] ; then
            args+=" --diamond"
        fi
        wrapper_phage_contigs_sorter_iPlant.pl \${args}
        """
}

workflow {
    qry = channel.fromPath(params.test_virsorter_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
    //qry.view { i -> "$i" }
    out = virsorter(qry)
    //out.view { i -> "$i" }
}

