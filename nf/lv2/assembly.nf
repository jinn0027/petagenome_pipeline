#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.assembly_filter_and_rename_memory = params.memory
params.assembly_filter_and_rename_threads = params.threads

params.assembly_get_length_memory = params.memory
params.assembly_get_length_threads = params.threads

params.assembly_get_stats_memory = params.memory
params.assembly_get_stats_threads = params.threads

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"
include { spades_assembler } from "${params.petagenomeDir}/nf/lv1/spades"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"

process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.assembly_filter_and_rename_memory}"
    def threads = "${params.assembly_filter_and_rename_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(id), path(read, arity: '1'), val(l_thre)
    output:
        tuple val(id), path("${id}/contig.${l_thre}.fa", arity: '0..*'), path("${id}/contig.name.txt")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${id}
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min ${l_thre} --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.${l_thre}.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 1000 --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.1000.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 5000 ${id}/contig.1000.fa > ${id}/contig.5000.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 10000 ${id}/contig.1000.fa > ${id}/contig.10000.fa
        """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.assembly_get_length_memory}"
    def threads = "${params.assembly_get_length_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(id), path(reads, arity: '0..*')
    output:
        tuple val(id), path("${id}/*.length.txt", arity: '0..*')
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
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.assembly_get_stats_memory}"
    def threads = "${params.assembly_get_stats_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
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
              #R --vanilla --slave --args \${i} 2 <  ${params.petagenomeDir}/scripts/R/stats.assembly.R > ${id}/\${outname}
              Rscript ${params.petagenomeDir}/scripts/R/stats.assembly.R \${i} 2 > ${id}/\${outname}
          else
              touch ${id}/\${outname}
          fi
        done
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// アセンブリ〜フィルタリング〜アノテーションDB作成の一連処理
workflow ASSEMBLY_SUB {
    take:
    p
    reads
    l_thre

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    // 1. SPAdes アセンブラの実行
    asm_raw = spades_assembler(in_ch)

    // scaffolds が存在する場合はそれを使い、無ければ contigs を選択
    asm = asm_raw.map { id, scaffolds, contigs ->
        tuple(id, (scaffolds && 0 < scaffolds.size()) ? scaffolds : contigs)
    }

    // 2. 配列長のフィルタリングとリネーム
    flt_in = p.combine(
        asm.map { id, contigs -> tuple(id, contigs, l_thre) }
    ).map { p_val, id, contigs, length_threshold ->
        tuple(p_val, id, contigs, length_threshold)
    }

    flt_raw = filter_and_rename(flt_in)

    // 後続プロセス用の入力作成 [ p_val, id, contigs ]
    flt_seqs = p.combine(
        flt_raw.map { id, contigs, name -> tuple(id, contigs) }
    ).map { p_val, id, contigs ->
        tuple(p_val, id, contigs)
    }

    // 3. 配列長計算、統計取得、BLAST DB作成
    len = get_length(flt_seqs)

    sts_in = p.combine(len).map { p_val, id, len_file ->
        tuple(p_val, id, len_file)
    }
    sts = get_stats(sts_in)

    blstdb = blast_makerefdb(flt_seqs)

    // 4. フィルタリング後のコンティグ配列を1ファイルずつ展開 (flatMap)
    flt = flt_raw.flatMap { id, contigs, name ->
        contigs.collect { c ->
            if (c && c.size() != 0) {
                return [c.getBaseName(), c]
            }
        }.findAll { it != null }
    }

    emit:
    asm    = asm
    flt    = flt
    len    = len
    sts    = sts
    blstdb = blstdb
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow ASSEMBLY_ALL {
    p      = createNullParamsChannel()
    reads  = createPairsChannel(params.test_assembly_reads)
    l_thre = params.test_assembly_l_thre

    out_ch = ASSEMBLY_SUB(p, reads, l_thre)

    // 各出力チャンネルの確認
    out_ch.asm.view    { i -> "ASM: $i" }
    out_ch.flt.view    { i -> "FLT: $i" }
    out_ch.len.view    { i -> "LEN: $i" }
    out_ch.sts.view    { i -> "STS: $i" }
    out_ch.blstdb.view { i -> "BLSTDB: $i" }
}

// デフォルトエントリーポイント
workflow {
    ASSEMBLY_ALL()
}
