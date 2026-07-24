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
params.mmseqs2_search_sort_results = 0 // Sort results: 0: no sorting, 1: sort by E-value (Alignment) or seq.id. (Hamming) [0]

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
params.mmseqs2_cluster_sort_results = 0 // Sort results: 0: no sorting, 1: sort by E-value (Alignment) or seq.id. (Hamming) [0]

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
params.mmseqs2_linclust_sort_results = 0 // Sort results: 0: no sorting, 1: sort by E-value (Alignment) or seq.id. (Hamming) [0]

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
        tuple val(p), val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${ref_id}
        mmseqs createdb \\
            --dbtype ${getParam(p, 'mmseqs2_ref_type')} \\
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
        tuple val(p), val(qry_id), path(qry, arity: '1')
    output:
        tuple val(qry_id), path("${qry_id}")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${qry_id}
        mmseqs createdb \\
            --dbtype ${getParam(p, 'mmseqs2_qry_type')} \\
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
        tuple val(p), val(ref_id), path(ref_db)
    output:
        tuple val(ref_id), path("${ref_id}/out.fasta"), path("${ref_id}/out.tsv")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${ref_id} tmp
        if [ "${getParam(p, 'mmseqs2_cluster_mode')}" == "cluster" ] ; then
            args="\\
                 -s ${getParam(p, 'mmseqs2_cluster_s')} \\
                 -k ${getParam(p, 'mmseqs2_cluster_k')} \\
                 -e ${getParam(p, 'mmseqs2_cluster_e')} \\
                 -c ${getParam(p, 'mmseqs2_cluster_c')} \\
                 --cov-mode ${getParam(p, 'mmseqs2_cluster_cov_mode')} \\
                 --min-seq-id ${getParam(p, 'mmseqs2_cluster_min_seq_id')} \\
                 --min-aln-len ${getParam(p, 'mmseqs2_cluster_min_aln_len')} \\
                 --max-seqs ${getParam(p, 'mmseqs2_cluster_max_seqs')} \\
                 --split ${getParam(p, 'mmseqs2_cluster_split')} \\
                 --split-mode ${getParam(p, 'mmseqs2_cluster_split_mode')} \\
                 --split-memory-limit ${getParam(p, 'mmseqs2_cluster_split_memory_limit')} \\
                 --sort-results ${getParam(p, 'mmseqs2_cluster_sort_results')} \\
                 "
        else
            args="\\
                 -k ${getParam(p, 'mmseqs2_linclust_k')} \\
                 -e ${getParam(p, 'mmseqs2_linclust_e')} \\
                 -c ${getParam(p, 'mmseqs2_linclust_c')} \\
                 --cov-mode ${getParam(p, 'mmseqs2_linclust_cov_mode')} \\
                 --min-seq-id ${getParam(p, 'mmseqs2_linclust_min_seq_id')} \\
                 --min-aln-len ${getParam(p, 'mmseqs2_linclust_min_aln_len')} \\
                 --split-memory-limit ${getParam(p, 'mmseqs2_linclust_split_memory_limit')} \\
                 --sort-results ${getParam(p, 'mmseqs2_linclust_sort_results')} \\
                 "
        fi
        mmseqs ${getParam(p, 'mmseqs2_cluster_mode')} \\
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
        tuple val(p), val(ref_id), path(ref_db), val(qry_id), path(qry_db)
    output:
        tuple val(ref_id), path("${qry_id}/out.*")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${qry_id} tmp
        args="\\
             --search-type ${getParam(p, 'mmseqs2_search_type')} \\
             -s ${getParam(p, 'mmseqs2_cluster_s')} \\
             -k ${getParam(p, 'mmseqs2_cluster_k')} \\
             -e ${getParam(p, 'mmseqs2_cluster_e')} \\
             -c ${getParam(p, 'mmseqs2_cluster_c')} \\
             --cov-mode ${getParam(p, 'mmseqs2_cluster_cov_mode')} \\
             --min-seq-id ${getParam(p, 'mmseqs2_cluster_min_seq_id')} \\
             --min-aln-len ${getParam(p, 'mmseqs2_cluster_min_aln_len')} \\
             --max-seqs ${getParam(p, 'mmseqs2_cluster_max_seqs')} \\
             --split ${getParam(p, 'mmseqs2_cluster_split')} \\
             --split-mode ${getParam(p, 'mmseqs2_cluster_split_mode')} \\
             --split-memory-limit ${getParam(p, 'mmseqs2_cluster_split_memory_limit')} \\
             --sort-results ${getParam(p, 'mmseqs2_search_sort_results')} \\
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

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// リファレンス DB 作成処理の本体
workflow BUILD_REF_DB_SUB {
    take:
    p
    ref

    main:
    db_input = p.combine(ref).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }
    ref_db = mmseqs2_makerefdb(db_input)

    emit:
    ref_db = ref_db
}

// クエリ DB 作成処理の本体
workflow BUILD_QRY_DB_SUB {
    take:
    p
    qry

    main:
    db_input = p.combine(qry).map { p_val, qry_id, qry_path ->
        tuple(p_val, qry_id, qry_path)
    }
    qry_db = mmseqs2_makeqrydb(db_input)

    emit:
    qry_db = qry_db
}

// MMseqs2 Search 処理の本体
workflow SEARCH_SUB {
    take:
    p
    ref_db
    qry_db

    main:
    // ref_db と qry_db を結合してフラット化
    in_ch = ref_db.combine(qry_db).map { ref_id, ref_path, qry_id, qry_path ->
        tuple(ref_id, ref_path, qry_id, qry_path)
    }

    // パラメータ p と結合してフラット化
    in_ch = p.combine(in_ch).map { p_val, ref_id, ref_path, qry_id, qry_path ->
        tuple(p_val, ref_id, ref_path, qry_id, qry_path)
    }

    out = mmseqs2_search(in_ch)

    emit:
    out = out
}

// MMseqs2 Cluster 処理の本体
workflow CLUSTER_SUB {
    take:
    p
    ref_db

    main:
    in_ch = p.combine(ref_db).map { p_val, ref_id, ref_path ->
        tuple(p_val, ref_id, ref_path)
    }

    out = mmseqs2_cluster(in_ch)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. DB 作成のみ実行 (-entry BUILD_DB_ONLY)
workflow BUILD_DB_ONLY {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.test_mmseqs2_ref)

    BUILD_REF_DB_SUB(p, ref)
}

// B. 作成済み DB を使って Search のみ実行 (-entry SEARCH_ONLY)
workflow SEARCH_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.test_mmseqs2_db)
    qry    = createSeqsChannel(params.test_mmseqs2_qry)

    qry_db = BUILD_QRY_DB_SUB(p, qry)
    out_ch = SEARCH_SUB(p, ref_db, qry_db.qry_db)

    out_ch.out.view { i -> "$i" }
}

// C. 作成済み DB を使って Cluster のみ実行 (-entry CLUSTER_ONLY)
workflow CLUSTER_ONLY {
    p      = createNullParamsChannel()
    ref_db = createSeqsChannel(params.test_mmseqs2_db)

    out_ch = CLUSTER_SUB(p, ref_db)
    out_ch.out.view { i -> "$i" }
}

// D. 条件分岐または一括実行 (デフォルト または -entry MMSEQS2_ALL)
workflow MMSEQS2_ALL {
    p   = createNullParamsChannel()
    ref = createSeqsChannel(params.test_mmseqs2_ref)

    ref_db_ch = BUILD_REF_DB_SUB(p, ref)

    if (params.test_mmseqs2_module == "search") {
        qry       = createSeqsChannel(params.test_mmseqs2_qry)
        qry_db_ch = BUILD_QRY_DB_SUB(p, qry)

        out_ch = SEARCH_SUB(p, ref_db_ch.ref_db, qry_db_ch.qry_db)
        out_ch.out.view { i -> "$i" }

    } else if (params.test_mmseqs2_module == "cluster") {
        out_ch = CLUSTER_SUB(p, ref_db_ch.ref_db)
        out_ch.out.view { i -> "$i" }
    }
}

// デフォルトエントリーポイント
workflow {
    MMSEQS2_ALL()
}
