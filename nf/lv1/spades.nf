#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. SPAdes 固有の上限値定義
def SPADES_EC_MAX_MEMORY         = 128 // GB (k-mer エラー訂正テーブル保持用)
def SPADES_EC_MAX_THREADS        = 16  // エラー訂正・Gzip 圧縮の並列飽和点

def SPADES_ASSEMBLY_MAX_MEMORY   = 250 // GB (De Bruijn Graph 構築用の高メモリ上限)
def SPADES_ASSEMBLY_MAX_THREADS  = 16  // グラフ走査・並列化効率の上限

// 3. 上限値による動的クリッピング
params.spades_spades_error_correction_memory             = Math.min(params.memory as Integer, SPADES_EC_MAX_MEMORY)
params.spades_spades_error_correction_threads            = Math.min(params.threads as Integer, SPADES_EC_MAX_THREADS)

params.spades_spades_error_correction_gzip_output_memory = Math.min(params.memory as Integer, SPADES_EC_MAX_MEMORY)
params.spades_spades_error_correction_gzip_output_threads= Math.min(params.threads as Integer, SPADES_EC_MAX_THREADS)

params.spades_spades_assembler_memory = Math.min(params.memory as Integer, SPADES_ASSEMBLY_MAX_MEMORY)
params.spades_spades_assembler_threads = Math.min(params.threads as Integer, SPADES_ASSEMBLY_MAX_THREADS)

params.spades_spades_e2e_memory = Math.min(params.memory as Integer, SPADES_ASSEMBLY_MAX_MEMORY)
params.spades_spades_e2e_threads = Math.min(params.threads as Integer, SPADES_ASSEMBLY_MAX_THREADS)

params.spades_e2e = false

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process spades_error_correction {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.spades_spades_error_correction_memory}"
    def threads = "${params.spades_spades_error_correction_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/corrected/paired/*.cor.fastq", arity: '0..2'), \
              path("${pair_id}/corrected/unpaired/*.cor.fastq", arity: '0..*')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir ${pair_id}
        spades.py \\
            --memory ${gb} \\
            --threads ${threads} \\
            --only-error-correction \\
            --disable-gzip-output \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        mkdir -p ${pair_id}/corrected/paired ${pair_id}/corrected/unpaired 
        mv ${pair_id}/corrected/*unpaired*.cor.fastq ${pair_id}/corrected/unpaired
        mv ${pair_id}/corrected/*.cor.fastq ${pair_id}/corrected/paired
        """
}

process spades_error_correction_gzip_output {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.spades_spades_error_correction_gzip_output_memory}"
    def threads = "${params.spades_spades_error_correction_gzip_output_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/corrected/paired/*.cor.fastq.gz", arity: '0..2'), \
              path("${pair_id}/corrected/unpaired/*.cor.fastq.gz", arity: '0..*')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir ${pair_id}
        spades.py \\
            --memory ${gb} \\
            --threads ${threads} \\
            --only-error-correction \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        mkdir -p ${pair_id}/corrected/paired ${pair_id}/corrected/unpaired 
        mv ${pair_id}/corrected/*unpaired*.cor.fastq.gz ${pair_id}/corrected/unpaired
        mv ${pair_id}/corrected/*.cor.fastq.gz ${pair_id}/corrected/paired
        """
}

process spades_assembler {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.spades_spades_assembler_memory}"
    def threads = "${params.spades_spades_assembler_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/scaffolds.fasta", arity: '0..*'), \
              path("${pair_id}/contigs.fasta", arity: '0..*')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        touch ${pair_id}/scaffolds.fasta ${pair_id}/contigs.fasta
        spades.py \\
            --memory ${gb} \\
            --threads ${threads} \\
            --only-assembler \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        """
}

process spades_e2e {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/spades/spades.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink'
    def gb = "${params.spades_spades_e2e_memory}"
    def threads = "${params.spades_spades_e2e_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), \
              path("${pair_id}/scaffolds.fasta", arity: '0..*'), \
              path("${pair_id}/contigs.fasta", arity: '0..*')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}
        spades.py \\
            --memory ${params.spades_spades_e2e_memory} \\
            --threads ${threads} \\
            --meta \\
            --pe1-1 ${reads[0]} \\
            --pe1-2 ${reads[1]} \\
            -o ${pair_id}
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// SPAdes End-to-End 処理の本体
workflow SPADES_E2E_SUB {
    take:
    p
    reads

    main:
    // [ p_val, pair_id, reads_path ] にフラット化
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = spades_e2e(in_ch)

    emit:
    out = out
}

// エラー修正 ＋ アセンブリ分割処理の本体
workflow SPADES_STEPWISE_SUB {
    take:
    p
    reads

    main:
    // エラー修正用入力のフラット化
    ec_in = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    // エラー修正実行と出力の整列 [ pair_id, paired ]
    ec_raw = spades_error_correction_gzip_output(ec_in)
    ec = ec_raw.map { pair_id, paired, unpaired ->
        tuple(pair_id, paired)
    }

    // アセンブリ用入力の結合とフラット化
    asm_in = p.combine(ec).map { p_val, pair_id, paired_path ->
        tuple(p_val, pair_id, paired_path)
    }

    out = spades_assembler(asm_in)

    emit:
    out = out
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. End-to-End 実行のみ (-entry SPADES_E2E_ONLY)
workflow SPADES_E2E_ONLY {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.spades_reads)

    out_ch = SPADES_E2E_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// B. エラー修正〜アセンブリ分割実行のみ (-entry SPADES_STEPWISE_ONLY)
workflow SPADES_STEPWISE_ONLY {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.spades_reads)

    out_ch = SPADES_STEPWISE_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

// C. 条件分岐または一括実行 (デフォルト または -entry SPADES_ALL)
workflow SPADES_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.spades_reads)

    if (params.spades_e2e) {
        out_ch = SPADES_E2E_SUB(p, reads)
        out_ch.out.view { i -> "$i" }
    } else {
        out_ch = SPADES_STEPWISE_SUB(p, reads)
        out_ch.out.view { i -> "$i" }
    }
}

// デフォルトエントリーポイント
workflow {
    SPADES_ALL()
}
