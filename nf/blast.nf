#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.blast_dbtype = "nucl"
params.blast_task = "megablast"
params.blast_num_alignments = "1"
params.blast_perc_identity = "95"
params.blast_evalue = "1e-20"
params.blast_outfmt = 6

process blast_makedb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/blast/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("db")
    script:
        """
        ref_=${ref}
        echo ${ref} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            ref_=\${ref_%%.gz}
            unpigz -c ${ref} > \${ref_}
        fi
        mkdir -p db
        makeblastdb \\
            -in \${ref_} \\
            -out db/${ref_id} \\
            -dbtype ${params.blast_dbtype} \\
            -parse_seqids
        """
}

process blastn {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/blast/${ref_id}/${qry_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(db, arity: '1'), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("out/*.tsv", arity: '1')
    script:
        """
        mkdir -p out
        blastn \\
            -task ${params.blast_task} \\
            -num_threads ${params.threads} \\
            -query ${qry} \\
            -db ${db}/${ref_id} \\
            -perc_identity ${params.blast_perc_identity} \\
            -evalue ${params.blast_evalue} \\
            -outfmt ${params.blast_outfmt} \\
            -num_alignments ${params.blast_num_alignments} \\
            -out out/${ref_id}_@_${qry_id}_out.tsv
        """
}

workflow {
    ref = channel.fromPath(params.test_blast_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }

    qry = channel.fromPath(params.test_blast_qry, checkIfExists: true)
        .map { qry_path -> tuple(qry_path.simpleName, qry_path) }
   
    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    db = blast_makedb(ref)
    //db.view { i -> "$i" }
    in = db.combine(qry)
    //in.view { i -> "$i" }
    out = blastn(in)
    out.view { i -> "$i" }
}

