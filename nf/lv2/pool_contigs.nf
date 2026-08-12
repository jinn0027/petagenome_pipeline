#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================
// 0. グローバルフォールバックと上限値（Clipping）定義
// ==========================================
params.memory  = params.memory  ?: 16
params.threads = params.threads ?: 4

// クラスタリングプロセスのデフォルト設定
params.pool_contigs_clustering_process = params.pool_contigs_clustering_process ?: "mmseqs2"

// プロセス固有の上限設定 (Clipping)
def MERGE_MAX_MEMORY          = 16
def MERGE_MAX_THREADS         = 8
def FILTER_RENAME_MAX_MEMORY  = 8
def FILTER_RENAME_MAX_THREADS = 4
def SUMMARIZE_MAX_MEMORY      = 8
def SUMMARIZE_MAX_THREADS     = 4
def GET_LEN_MAX_MEMORY        = 8
def GET_LEN_MAX_THREADS       = 4
def GET_STATS_MAX_MEMORY      = 8
def GET_STATS_MAX_THREADS     = 2

// リソース割り当てと上限適応
params.pool_contigs_mergs_contigs_memory    = Math.min((params.pool_contigs_mergs_contigs_memory    ?: params.memory)  as Integer, MERGE_MAX_MEMORY)
params.pool_contigs_mergs_contigs_threads   = Math.min((params.pool_contigs_mergs_contigs_threads   ?: params.threads) as Integer, MERGE_MAX_THREADS)

params.pool_contigs_filter_and_rename_memory  = Math.min((params.pool_contigs_filter_and_rename_memory  ?: params.memory)  as Integer, FILTER_RENAME_MAX_MEMORY)
params.pool_contigs_filter_and_rename_threads = Math.min((params.pool_contigs_filter_and_rename_threads ?: params.threads) as Integer, FILTER_RENAME_MAX_THREADS)

params.pool_contigs_summarize_name_memory     = Math.min((params.pool_contigs_summarize_name_memory     ?: params.memory)  as Integer, SUMMARIZE_MAX_MEMORY)
params.pool_contigs_summarize_name_threads    = Math.min((params.pool_contigs_summarize_name_threads    ?: params.threads) as Integer, SUMMARIZE_MAX_THREADS)

params.pool_contigs_get_length_memory       = Math.min((params.pool_contigs_get_length_memory       ?: params.memory)  as Integer, GET_LEN_MAX_MEMORY)
params.pool_contigs_get_length_threads      = Math.min((params.pool_contigs_get_length_threads      ?: params.threads) as Integer, GET_LEN_MAX_THREADS)

params.pool_contigs_get_stats_memory        = Math.min((params.pool_contigs_get_stats_memory        ?: params.memory)  as Integer, GET_STATS_MAX_MEMORY)
params.pool_contigs_get_stats_threads       = Math.min((params.pool_contigs_get_stats_threads       ?: params.threads) as Integer, GET_STATS_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"
include { cdhit_est } from "${params.petagenomeDir}/nf/lv1/cdhit"
include { mmseqs2_makerefdb; mmseqs2_cluster } from "${params.petagenomeDir}/nf/lv1/mmseqs2"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"

// ==========================================
// 1. プロセス定義
// ==========================================

process merge_contigs {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.pool_contigs_mergs_contigs_memory}"
    def threads = "${params.pool_contigs_mergs_contigs_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(contigs, arity: '1..*', stageAs: "?/*")

    output:
        tuple val(id), path("${id}/merged_contig.fa"), path("${id}/contig_list.txt")

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    touch ${id}/merged_contig.fa
    touch ${id}/contig_list.txt
    contigs_=( ${contigs} )
    for contig in \${contigs_[@]}
    do
        contig_=\${contig}
        echo \${contig} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            contig_=\${contig_%%.gz}
            unpigz -c \${contig} > \${contig_}
        fi
        cat \${contig_} >> ${id}/merged_contig.fa
        echo \${contig_} >> ${id}/contig_list.txt
    done
    """
}

process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.pool_contigs_filter_and_rename_memory}"
    def threads = "${params.pool_contigs_filter_and_rename_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(read, arity: '1'), val(l_thre)

    output:
        tuple val(id), path("${id}/contig.${l_thre}.fa"), path("${id}/contig.name.txt")

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
         --min ${l_thre} --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.${l_thre}.fa
    """
}

process summarize_name {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.pool_contigs_summarize_name_memory}"
    def threads = "${params.pool_contigs_summarize_name_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(name, arity: '1')
        tuple val(id), path(clstr, arity: '1')

    output:
        tuple val(id), path("${id}/${id}.name.txt"), path("${id}/*")

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    if [ "${getParam(p, 'pool_contigs_clustering_process')}" = "mmseqs2" ] ; then
        cp -f ${clstr} ${id}/${id}.name_
    else
        ruby ${params.petagenomeDir}/scripts/Ruby/parse.cdhit_clstr.rb -i ${clstr} --include_rep > ${id}/${id}.name_
    fi
    ruby ${params.petagenomeDir}/scripts/Ruby/join_with_tab.rb ${id}/${id}.name_ 2 ${name} 2 | \
         awk -F '\\t' '{OFS="\\t"} {print \$2,\$1,\$3}' > ${id}/${id}.name.txt
    awk '{print(\$1)}' ${id}/${id}.name.txt | sort | uniq > ${id}/${id}.samples.txt
    while read sample
    do
        awk -F "\\t" -v sample=\${sample} \
            '{OFS="\\t"} { if (\$1 ~ sample) print \$0 }' \
            ${id}/${id}.name.txt \
            > ${id}/${id}.\${sample}.name.txt
    done < ${id}/${id}.samples.txt
    rm -f ${id}/${id}.name_
    """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.pool_contigs_get_length_memory}"
    def threads = "${params.pool_contigs_get_length_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(reads, arity: '1..*')

    output:
        tuple val(id), path("${id}/*.length.txt")

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    reads_=( ${reads} )
    for i in \${reads_[@]}
    do
        awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' \${i} | \
        python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fasta \
        > ${id}/\$(basename \$i).length.txt
    done
    """
}

process get_stats {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.pool_contigs_get_stats_memory}"
    def threads = "${params.pool_contigs_get_stats_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(lengths, arity: '1..*')

    output:
        tuple val(id), path("${id}/*.stats.txt")

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    lengths_=( ${lengths} )
    for i in \${lengths_[@]}
    do
        outname=\$(basename \${i} | sed "s#length#stats#")
        n=\$(cat \${i} | wc -l)
        if [ \${n} -gt 0 ] ; then
            Rscript ${params.petagenomeDir}/scripts/R/stats.assembly.R \${i} 2 > ${id}/\${outname}
        else
            touch ${id}/\${outname}
        fi
    done
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow POOL_CONTIGS_SUB {
    take:
    p
    contigs
    l_thre

    main:
    // [ p_val, seq_id, seq_path ] に整形
    in_ch = p.combine(contigs).map { p_val, seq_id, seq_path ->
        tuple(p_val, seq_id, seq_path)
    }

    // 1. アセンブリの結合
    merged = merge_contigs(in_ch)

    // 結合済み FASTA の抽出 [ p_val, id, fasta ]
    merged_fasta = p.combine(
        merged.map { id, fasta, list -> tuple(id, fasta) }
    ).map { p_val, id, fasta ->
        tuple(p_val, id, fasta)
    }

    // 2. クラスタリングによる重複除去 (MMseqs2 または CD-HIT)
    if (params.pool_contigs_clustering_process == "mmseqs2") {
        merged_db = mmseqs2_makerefdb(merged_fasta)
        clust_in = p.combine(merged_db).map { p_val, id, db_path ->
            tuple(p_val, id, db_path)
        }
        clust = mmseqs2_cluster(clust_in)
    } else {
        clust = cdhit_est(merged_fasta)
    }

    // 3. リネームと配列長フィルタリング (L >= l_thre)
    flt_in = p.combine(
        clust.map { id, fasta, clstr -> tuple(id, fasta, l_thre) }
    ).map { p_val, id, fasta, length_threshold ->
        tuple(p_val, id, fasta, length_threshold)
    }
    flt = filter_and_rename(flt_in)

    // フィルタリング後の配列の抽出 [ p_val, id, fasta ]
    flt_seqs = p.combine(
        flt.map { id, fasta, name -> tuple(id, fasta) }
    ).map { p_val, id, fasta ->
        tuple(p_val, id, fasta)
    }

    // 4. 名前の要約テーブル作成
    flt_names = p.combine(
        flt.map { id, fasta, name -> tuple(id, name) }
    ).map { p_val, id, name ->
        tuple(p_val, id, name)
    }
    clust_clstr = clust.map { id, fasta, clstr -> tuple(id, clstr) }

    name = summarize_name(flt_names, clust_clstr)

    // 5. 配列長計算 ＋ 統計取得 ＋ BLAST DB作成
    len = get_length(flt_seqs)

    sts_in = p.combine(len).map { p_val, id, len_file ->
        tuple(p_val, id, len_file)
    }
    sts = get_stats(sts_in)

    blstdb = blast_makerefdb(flt_seqs)

    emit:
    merged = merged
    clust  = clust
    flt    = flt
    name   = name
    len    = len
    sts    = sts
    blstdb = blstdb
}

// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow POOL_CONTIGS_ALL {
    p        = createNullParamsChannel()
    contigs  = createSeqsChannel(params.pool_contigs_contigs)
    l_thresh = params.pool_contigs_l_thre

    out_ch = POOL_CONTIGS_SUB(p, contigs, l_thresh)

    out_ch.merged.view { i -> "MERGED: $i" }
    out_ch.clust.view  { i -> "CLUST: $i" }
    out_ch.flt.view    { i -> "FLT: $i" }
    out_ch.name.view   { i -> "NAME: $i" }
    out_ch.len.view    { i -> "LEN: $i" }
    out_ch.sts.view    { i -> "STS: $i" }
    out_ch.blstdb.view { i -> "BLSTDB: $i" }
}

workflow {
    POOL_CONTIGS_ALL()
}