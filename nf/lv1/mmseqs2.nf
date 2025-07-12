#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.mmseqs2_mmseqs2_makerefdb_memory = params.memory
params.mmseqs2_mmseqs2_makerefdb_threads = params.threads

params.mmseqs2_mmseqs2_makeqrydb_memory = params.memory
params.mmseqs2_mmseqs2_makeqrydb_threads = params.threads

params.mmseqs2_mmseqs2_cluster_memory = params.memory
params.mmseqs2_mmseqs2_cluster_threads = params.threads

params.mmseqs2_mmseqs2_search_memory = params.memory
params.mmseqs2_mmseqs2_search_threads = params.threads

params.mmseqs2_ref_type = "0" //Database type 0: auto, 1: amino acid 2: nucleotides [0]
params.mmseqs2_qry_type = "0" //Database type 0: auto, 1: amino acid 2: nucleotides [0]
params.mmseqs2_search_type = "0" // Search type 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment [0]
params.mmseqs2_cluster_mode = "cluster" // cluster or linclust

process mmseqs2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.mmseqs2_mmseqs2_makerefdb_memory} GB"
    cpus "${params.mmseqs2_mmseqs2_makerefdb_threads}"

    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        mkdir -p ${ref_id}
        mmseqs createdb \\
            --dbtype ${params.mmseqs2_ref_type} \\
            ${ref} \\
            ${ref_id}/ref
        """
}

process mmseqs2_makeqrydb {
    tag "${qry_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.mmseqs2_mmseqs2_makeqrydb_memory} GB"
    cpus "${params.mmseqs2_mmseqs2_makeqrydb_threads}"

    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id), path("${qry_id}")
    script:
        """
        mkdir -p ${qry_id}
        mmseqs createdb \\
            --dbtype ${params.mmseqs2_qry_type} \\
            ${qry} \\
            ${qry_id}/qry
        """
}

process mmseqs2_cluster {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.mmseqs2_mmseqs2_cluster_memory} GB"
    cpus "${params.mmseqs2_mmseqs2_cluster_threads}"

    input:
        tuple val(ref_id), path(ref_db)
    output:
        tuple val(ref_id), path("${ref_id}/out.fasta"), path("${ref_id}/out.tsv")
    script:
        """
        mkdir -p ${ref_id} tmp
        mmseqs ${params.mmseqs2_cluster_mode} \\
            --threads ${params.mmseqs2_mmseqs2_cluster_threads} \\
            ${ref_db}/ref \\
            ${ref_id}/clu \\
            tmp
        mmseqs createtsv \\
            --threads ${params.threads} \\
            ${ref_db}/ref \\
            ${ref_db}/ref \\
            ${ref_id}/clu \\
            ${ref_id}/out.tsv
        mmseqs result2repseq \\
            --threads ${params.threads} \\
            ${ref_db}/ref \\
            ${ref_id}/clu \\
            ${ref_id}/clu_rep
        mmseqs result2flat \\
            ${ref_db}/ref \\
            ${ref_db}/ref \\
            ${ref_id}/clu_rep \\
            ${ref_id}/out.fasta \\
            --use-fasta-header
        """
}

process mmseqs2_search {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory "${params.mmseqs2_mmseqs2_search_memory} GB"
    cpus "${params.mmseqs2_mmseqs2_search_threads}"

    input:
        tuple val(ref_id), path(ref_db), val(qry_id), path(qry_db)
    output:
        tuple val(ref_id), path("${qry_id}/out.*")
    script:
        """
        mkdir -p ${qry_id} tmp
        mmseqs search \\
            --threads ${params.mmseqs2_mmseqs2_search_threads} \\
            --search-type ${params.mmseqs2_search_type} \\
            ${qry_db}/qry \\
            ${ref_db}/ref \\
            ${qry_id}/out \\
            tmp
        """
}

workflow {
    ref = channel.fromPath(params.test_mmseqs2_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }
    //ref.view { i -> "$i" }

    ref_db = mmseqs2_makerefdb(ref)
    //ref_db.view { i -> "$i" }

    if (params.test_mmseqs2_module == "search") {

        qry = channel.fromPath(params.test_mmseqs2_qry, checkIfExists: true)
            .map { qry_path -> tuple(qry_path.simpleName, qry_path) }

        //qry.view { i -> "$i" }

        qry_db = mmseqs2_makeqrydb(qry)
        //qry_db.view { i -> "$i" }

        in = ref_db.combine(qry_db)
        //in.view { i -> "$i" }

        out = mmseqs2_search(in)
        //out.view { i -> "$i" }

    } else if (params.test_mmseqs2_module == "cluster") {

        out = mmseqs2_cluster(ref_db)
        //out.view { i -> "$i" }

    }
}

