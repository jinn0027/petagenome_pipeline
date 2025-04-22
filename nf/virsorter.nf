#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.virsorter_data_dir = "/opt/VirSorter/virsorter-data"
params.virsorter_db = "refseq"
//params.virsorter_db = "virome"
params.virsorter_aligner = "blast"
//params.virsorter_aligner = "diamond"

process virsorter {
    tag "${qry_id}"
    container = "${params.petagenomeDir}/modules/virsorter/virsorter.sif"
    publishDir "${params.output}/virsorter/${params.virsorter_db}/${qry_id}", mode: 'copy'
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id) , path("${params.virsorter_aligner}/*.csv", arity: '1')
    script:
        """
        qry_=${qry}
        echo ${qry} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            qry_=\${qry_%%.gz}
            unpigz -c ${qry} > \${qry_}
        fi
        db=${params.virsorter_db}
        if [ "\${db}" = "refseq" ] ; then
            db="1"
        elif [ "\${db}" = "virome" ] ; then
            db="2"
        fi
        args="\\
            --db \${db} \\
            --wdir ${params.virsorter_aligner} \\
            --data-dir ${params.virsorter_data_dir} \\
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
    qry.view { i -> "$i" }
    virsorter = virsorter(qry)
    //virsorter.view { i -> "$i" }
}

