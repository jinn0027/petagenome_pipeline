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
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"
include { spades_assembler } from "${params.petagenomeDir}/nf/lv1/spades.nf"
include { MEGAHIT_SUB      } from "${params.petagenomeDir}/nf/lv1/megahit.nf"


// ==========================================
// 1. プロセス定義
// ==========================================

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
        tuple val(id), path(read, arity: '1'), val(l_thre)
    output:
        tuple val(id), path("${id}/contig.${l_thre}.fa", arity: '0..*'), path("${id}/contig.name.txt")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${id}
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min ${l_thre} --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.${l_thre}.fa
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
        tuple val(id), path(reads, arity: '0..*')
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
        tuple val(id), path(lengths, arity: '1..*')
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
// 2. サブワークフロー（再利用可能な処理の本体）
// ==========================================

workflow ASSEMBLY_SUB {
    take:
    p
    reads
    l_thre

    main:
    def tool = params.containsKey('assembler') ? params.assembler : 'megahit'

    // 1. アセンブラの分岐実行
    if (tool == 'spades') {
        in_ch = p.combine(reads).map { p_val, pair_id, reads_path -> tuple(p_val, pair_id, reads_path) }
        asm_raw = spades_assembler(in_ch)
        asm = asm_raw.map { id, scaffolds, contigs -> tuple(id, (scaffolds && 0 < scaffolds.size()) ? scaffolds : contigs) }
    } else {
        // MEGAHIT_SUB の出力は [pair_id, contigs]
        asm_raw = MEGAHIT_SUB(p, reads)
        asm = asm_raw.out
    }

    // 2. 配列長のフィルタリングとリネーム
    flt_in = asm.map { id, contigs -> 
        tuple(id, contigs, l_thre) 
    }
    flt_raw = filter_and_rename(flt_in)
    flt_seqs = flt_raw.map { id, contigs, name -> 
        tuple(id, contigs) 
    }

    // 3. 配列長計算、統計取得
    len = get_length(flt_seqs)
    sts = get_stats(len)

    emit:
    asm      = asm
    flt_seqs = flt_seqs
    len      = len
    sts      = sts
}


// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================
workflow ASSEMBLY_ALL {
    p      = createNullParamsChannel()
    reads  = createPairsChannel(params.assembly_reads)
    l_thre = params.containsKey('assembly_l_thre') ? params.assembly_l_thre : 1000

    out_ch = ASSEMBLY_SUB(p, reads, l_thre)

    out_ch.asm.view     { i -> "ASM: $i" }
    out_ch.flt_seq.view { i -> "FLT_SEQ: $i" }
    out_ch.len.view     { i -> "LEN: $i" }
    out_ch.sts.view     { i -> "STS: $i" }
}

workflow {
    ASSEMBLY_ALL()
}
