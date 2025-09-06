#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.blast_blast_makerefdb_memory = params.memory
params.blast_blast_makerefdb_threads = params.threads

params.blast_blastn_memory = params.memory
params.blast_blastn_threads = params.threads

params.blast_dbtype = "nucl"
params.blast_task = "megablast"
params.blast_num_alignments = "1"
params.blast_perc_identity = "95"
params.blast_evalue = "1e-20"
params.blast_outfmt = 6

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process blast_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.blast_blast_makerefdb_memory}"
    def threads = "${params.blast_blast_makerefdb_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref, arity: '1')
        val(p)
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}"
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
            -dbtype ${getParam(p, 'blast_dbtype')} \\
            -parse_seqids
        """
}

process blastn {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/blast/blast.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.blast_blastn_memory}"
    def threads = "${params.blast_blastn_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(qry_id), path(qry, arity: '1')
        val(p)
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/out.tsv", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        qry_=${qry}
        echo ${qry} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            qry_=\${qry_%%.gz}
            unpigz -c ${qry} > \${qry_}
        fi
        mkdir -p ${qry_id}
        blastn \\
            -task ${params.blast_task} \\
            -num_threads ${threads} \\
            -query \${qry_} \\
            -db ${ref_db}/ref \\
            -perc_identity ${getParam(p, 'blast_perc_identity')} \\
            -evalue ${getParam(p, 'blast_evalue')} \\
            -outfmt ${getParam(p, 'blast_outfmt')} \\
            -num_alignments ${getParam(p, 'blast_num_alignments')} \\
            -out ${qry_id}/out.tsv
        """
}

workflow {
    ref = createSeqsChannel(params.test_blast_ref)
    qry = createSeqsChannel(params.test_blast_qry)

    //ref.view { i -> "$i" }
    //qry.view { i -> "$i" }

    ref_db = blast_makerefdb(ref, createNullParamsChannel())
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(qry)
    //in.view { i -> "$i" }
    out = blastn(in, createNullParamsChannel())
    out.view { i -> "$i" }
}

