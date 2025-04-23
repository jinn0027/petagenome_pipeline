#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.virsorter_data_dir = "/opt/VirSorter/virsorter-data"
params.virsorter_db_type = "refseq"
//params.virsorter_db = "virome"
params.virsorter_aligner = "blast"
//params.virsorter_aligner = "diamond"

process virsorter {
    tag "${qry_id}"
    def local_db = "/opt/VirSorter/virsorter-data"
    def local_mga = "/opt/VirSorter/mga_linux_ia64"
    container = "${params.petagenomeDir}/modules/virsorter/virsorter.sif"
    containerOptions "-B ${db}:${local_db} -B ${mga}:${local_mga}"
    publishDir "${params.output}/virsorter/${params.virsorter_db}/${qry_id}", mode: 'copy'
    input:
        tuple path(db, arity: '1'), path(mga, arity: '1'), val(qry_id), path(qry, arity: '1')
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
        db_type=${params.virsorter_db_type}
        if [ "\${db_type}" = "refseq" ] ; then
            db_type="1"
        elif [ "\${db_type}" = "virome" ] ; then
            db_type="2"
        fi
        args="\\
            --db \${db_type} \\
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
    db = channel.fromPath(params.virsorter_db, checkIfExists: true)
    mga = channel.fromPath(params.virsorter_mga, checkIfExists: true)
    qry = channel.fromPath(params.test_virsorter_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
    in = db.combine(mga).combine(qry)
    //in.view { i -> "$i" }
    out = virsorter(in)
    //out.view { i -> "$i" }
}

