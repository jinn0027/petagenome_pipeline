#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.blast_blast_makerefdb__memory = params.memory
params.blast_blast_makerefdb_threads = params.threads

params.blast_blastn_memory = params.memory
params.blast_blastn_threads = params.threads

params.blast_dbtype = "nucl"
params.blast_task = "megablast"
params.blast_num_alignments = "1"
params.blast_perc_identity = "95"
params.blast_evalue = "1e-20"
params.blast_outfmt = 6

process blast_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory params.blast_blast_makerefdb__memory
    cpus params.blast_blast_makerefdb_threads

    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        ref_=${ref}
        echo ${ref} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            ref_=\${ref_%%.gz}
            unpigz -c ${ref} > \${ref_}
        fi
        mkdir -p ${ref_id}
        makeblastdb \\
            -in \${ref_} \\
            -out ${ref_id}/ref \\
            -dbtype ${params.blast_dbtype} \\
            -parse_seqids
        """
}

process blastn {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory params.blast_blastn_memory
    cpus params.blast_blastn_threads

    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.tsv", arity: '1')
    script:
        """
        mkdir -p ${qry_id}
        blastn \\
            -task ${params.blast_task} \\
            -num_threads ${params.blast_blastn_threads} \\
            -query ${qry} \\
            -db ${ref_db}/ref \\
            -perc_identity ${params.blast_perc_identity} \\
            -evalue ${params.blast_evalue} \\
            -outfmt ${params.blast_outfmt} \\
            -num_alignments ${params.blast_num_alignments} \\
            -out ${qry_id}/out.tsv
        """
}

workflow {
    ref = channel.fromPath(params.test_blast_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_blast_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = blast_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = blastn(in)
    out.view { i -> "$i" }
}

