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

params.mmseqs2_ref_type = 0 //Database type 0: auto, 1: amino acid 2: nucleotides [0]
params.mmseqs2_qry_type = 0 //Database type 0: auto, 1: amino acid 2: nucleotides [0]

//=== search params
params.mmseqs2_search_type = 0 // Search type 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment [0]
params.mmseqs2_search_s = 5.7 // Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [5.700]
params.mmseqs2_search_k = 15 // k-mer length (0: automatically set to optimum) [0]->[15]?
params.mmseqs2_search_e = "1.000e-03"
params.mmseqs2_search_c = 0.0 // List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]->[0.0]?
params.mmseqs2_search_cov_mode = 0 // 0: coverage of query and target
                                   // 1: coverage of target
                                   // 2: coverage of query
                                   // 3: target seq. length has to be at least x% of query length
                                   // 4: query seq. length has to be at least x% of target length
                                   // 5: short seq. needs to be at least x% of the other seq. length [0]
params.mmseqs2_search_min_seq_id = 0.0 // List matches above this sequence identity (for clustering) (range 0.0-1.0) [0.000]
params.mmseqs2_search_min_aln_len = 0 // Minimum alignment length (range 0-INT_MAX) [0]
params.mmseqs2_search_split = 0 // Split input into N equally distributed chunks. 0: set the best split automatically [0]
params.mmseqs2_search_split_mode = 2 // 0: split target db; 1: split query db; 2: auto, depending on main memory [2]
params.mmseqs2_search_split_memory_limit = 0 // Set max memory per split. E.g. 800B, 5K, 10M, 1G. Default (0) to all available system memory [0]
params.mmseqs2_search_max_seqs = 300 // Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [300]

//=== cluster params
params.mmseqs2_cluster_mode = "cluster" // cluster or linclust

//=== cluster [cluster]
params.mmseqs2_cluster_s = 4.0 // Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
params.mmseqs2_cluster_k = 15 // k-mer length (0: automatically set to optimum) [0]->[15]?
params.mmseqs2_cluster_e = "1.000e-03"
params.mmseqs2_cluster_c = 0.8 // List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]
params.mmseqs2_cluster_cov_mode = 0 // 0: coverage of query and target
                                    // 1: coverage of target
                                    // 2: coverage of query
                                    // 3: target seq. length has to be at least x% of query length
                                    // 4: query seq. length has to be at least x% of target length
                                    // 5: short seq. needs to be at least x% of the other seq. length [0]
params.mmseqs2_cluster_min_seq_id = 0.0 // List matches above this sequence identity (for clustering) (range 0.0-1.0) [0.000]
params.mmseqs2_cluster_min_aln_len = 0 // Minimum alignment length (range 0-INT_MAX) [0]
params.mmseqs2_cluster_split = 0 // Split input into N equally distributed chunks. 0: set the best split automatically [0]
params.mmseqs2_cluster_split_mode = 2 // 0: split target db; 1: split query db; 2: auto, depending on main memory [2]
params.mmseqs2_cluster_split_memory_limit = 0 // Set max memory per split. E.g. 800B, 5K, 10M, 1G. Default (0) to all available system memory [0]
params.mmseqs2_cluster_max_seqs = 20 // Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [20]

//=== cluster [linclust]
params.mmseqs2_linclust_k = 0 // k-mer length (0: automatically set to optimum) [0]
params.mmseqs2_linclust_e = "1.000e-03"
params.mmseqs2_linclust_c = 0.8 // List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]
params.mmseqs2_linclust_cov_mode = 0 // 0: coverage of query and target
                                     // 1: coverage of target
                                     // 2: coverage of query
                                     // 3: target seq. length has to be at least x% of query length
                                     // 4: query seq. length has to be at least x% of target length
                                     // 5: short seq. needs to be at least x% of the other seq. length [0]
params.mmseqs2_linclust_min_seq_id = 0.0 // List matches above this sequence identity (for clustering) (range 0.0-1.0) [0.900]
params.mmseqs2_linclust_min_aln_len = 0 // Minimum alignment length (range 0-INT_MAX) [0]
params.mmseqs2_linclust_split_memory_limit = 0 // Set max memory per split. E.g. 800B, 5K, 10M, 1G. Default (0) to all available system memory [0]

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"

process mmseqs2_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/mmseqs2/mmseqs2.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.mmseqs2_mmseqs2_makerefdb_memory}"
    def threads = "${params.mmseqs2_mmseqs2_makerefdb_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}"
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
    def gb = "${params.mmseqs2_mmseqs2_makeqrydb_memory}"
    def threads = "${params.mmseqs2_mmseqs2_makeqrydb_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id), path("${qry_id}")
    script:
        """
        echo "${processProfile(task)}"
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
    def gb = "${params.mmseqs2_mmseqs2_cluster_memory}"
    def threads = "${params.mmseqs2_mmseqs2_cluster_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref_db)
    output:
        tuple val(ref_id), path("${ref_id}/out.fasta"), path("${ref_id}/out.tsv")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${ref_id} tmp
        if [ "${params.mmseqs2_cluster_mode}" == "cluster" ] ; then
            args="\\
                 -s ${params.mmseqs2_cluster_s} \\
                 -k ${params.mmseqs2_cluster_k} \\
                 -e ${params.mmseqs2_cluster_e} \\
                 -c ${params.mmseqs2_cluster_c} \\
                 --cov-mode ${params.mmseqs2_cluster_cov_mode} \\
                 --min-seq-id ${params.mmseqs2_cluster_min_seq_id} \\
                 --min-aln-len ${params.mmseqs2_cluster_min_aln_len} \\
                 --max-seqs ${params.mmseqs2_cluster_max_seqs} \\
                 --split ${params.mmseqs2_cluster_split} \\
                 --split-mode ${params.mmseqs2_cluster_split_mode} \\
                 --split-memory-limit ${params.mmseqs2_cluster_split_memory_limit} \\
                 "
        else
            args="\\
                 -k ${params.mmseqs2_linclust_k} \\
                 -e ${params.mmseqs2_linclust_e} \\
                 -c ${params.mmseqs2_linclust_c} \\
                 --cov-mode ${params.mmseqs2_linclust_cov_mode} \\
                 --min-seq-id ${params.mmseqs2_linclust_min_seq_id} \\
                 --min-aln-len ${params.mmseqs2_linclust_min_aln_len} \\
                 --split-memory-limit ${params.mmseqs2_linclust_split_memory_limit} \\
                 "
        fi
        mmseqs ${params.mmseqs2_cluster_mode} \\
            --threads ${threads} \\
            \${args} \\
            ${ref_db}/ref \\
            ${ref_id}/clu \\
            tmp
        mmseqs createtsv \\
            --threads ${threads} \\
            ${ref_db}/ref \\
            ${ref_db}/ref \\
            ${ref_id}/clu \\
            ${ref_id}/out.tsv
        mmseqs result2repseq \\
            --threads ${threads} \\
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
    def gb = "${params.mmseqs2_mmseqs2_search_memory}"
    def threads = "${params.mmseqs2_mmseqs2_search_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(ref_id), path(ref_db), val(qry_id), path(qry_db)
    output:
        tuple val(ref_id), path("${qry_id}/out.*")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id} tmp
        args="\\
             --search-type ${params.mmseqs2_search_type} \\
             -s ${params.mmseqs2_cluster_s} \\
             -k ${params.mmseqs2_cluster_k} \\
             -e ${params.mmseqs2_cluster_e} \\
             -c ${params.mmseqs2_cluster_c} \\
             --cov-mode ${params.mmseqs2_cluster_cov_mode} \\
             --min-seq-id ${params.mmseqs2_cluster_min_seq_id} \\
             --min-aln-len ${params.mmseqs2_cluster_min_aln_len} \\
             --max-seqs ${params.mmseqs2_cluster_max_seqs} \\
             --split ${params.mmseqs2_cluster_split} \\
             --split-mode ${params.mmseqs2_cluster_split_mode} \\
             --split-memory-limit ${params.mmseqs2_cluster_split_memory_limit} \\
                 "
        mmseqs search \\
            --threads ${threads} \\
            \${args} \\
            ${qry_db}/qry \\
            ${ref_db}/ref \\
            ${qry_id}/out \\
            tmp
        """
}

workflow {
    p = createNullParamsChannel()
    ref = createSeqsChannel(params.test_mmseqs2_ref)
    //ref.view { i -> "$i" }

    ref_db = mmseqs2_makerefdb(ref)
    //ref_db.view { i -> "$i" }

    if (params.test_mmseqs2_module == "search") {
        qry = createSeqsChannel(params.test_mmseqs2_qry)
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

