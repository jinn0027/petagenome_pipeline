#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. PRINSEQ 固有の上限値定義
// ※ prinseq-lite.pl は Single Thread 動作のため MAX_THREADS は 1〜2 が推奨されます
def PRINSEQ_MAX_MEMORY  = 8 // GB
def PRINSEQ_MAX_THREADS = 1 // prinseq-lite.pl はマルチスレッド非対応のため 1 に設定

// 3. 上限値による動的クリッピング
params.prinseq_prinseq_memory  = Math.min(params.memory as Integer, PRINSEQ_MAX_MEMORY)
params.prinseq_prinseq_threads = Math.min(params.threads as Integer, PRINSEQ_MAX_THREADS)

params.prinseq_trim_right = 10
params.prinseq_trim_left = 10
params.prinseq_qual_right = 20
params.prinseq_qual_left = 20
params.prinseq_qual_window = 20
params.prinseq_min_len = 75
params.prinseq_derep = 1
params.prinseq_lc_method = "dust"
params.prinseq_lc_threshold = 7
params.prinseq_trim_ns_right = 1
params.prinseq_ns_max_n = 0

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process prinseq {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/prinseq/prinseq.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    cpus   { params.executor == "sge" ? null : params.prinseq_prinseq_threads }
    memory { params.executor == "sge" ? null : "${params.prinseq_prinseq_memory} GB" }
    clusterOptions "${clusterOptions(params.executor, params.prinseq_prinseq_memory, params.prinseq_prinseq_threads, label)}"

    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(pair_id), path("${pair_id}")
    script:
        def gb = params.prinseq_prinseq_memory
        def threads = params.prinseq_prinseq_threads
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        read0=${reads[0]}
        read1=${reads[1]}

        # gzip判定と解凍（unpigzがあれば並列解凍、なければgzipを使用）
        if [[ "\${read0}" =~ \\.gz\$ ]]; then
            read0=\${read0%%.gz}
            read1=\${read1%%.gz}
            
            DECOMPRESSOR="gzip -dc"
            if command -v unpigz >/dev/null 2>&1; then
                DECOMPRESSOR="unpigz -c -p ${task.cpus ?: 1}"
            fi

            \$DECOMPRESSOR ${reads[0]} > \${read0}
            \$DECOMPRESSOR ${reads[1]} > \${read1}
            IS_DECOMPRESSED=1
        fi

        prinseq-lite.pl \\
            -trim_left ${getParam(p, 'prinseq_trim_left')} \\
            -trim_right ${getParam(p, 'prinseq_trim_right')} \\
            -trim_qual_left ${getParam(p, 'prinseq_qual_left')} \\
            -trim_qual_right ${getParam(p, 'prinseq_qual_right')} \\
            -trim_qual_window ${getParam(p, 'prinseq_qual_window')} \\
            -min_len ${getParam(p, 'prinseq_min_len')} \\
            -derep ${getParam(p, 'prinseq_derep')} \\
            -lc_method ${getParam(p, 'prinseq_lc_method')} \\
            -lc_threshold ${getParam(p, 'prinseq_lc_threshold')} \\
            -trim_ns_right ${getParam(p, 'prinseq_trim_ns_right')} \\
            -ns_max_n ${getParam(p, 'prinseq_ns_max_n')} \\
            -out_good ${pair_id}/good \\
            -out_bad ${pair_id}/bad \\
            -fastq \${read0} \\
            -fastq2 \${read1}

        # 展開した一時FASTQファイルのクリーンアップ（ディスク空き容量の確保）
        if [ "\${IS_DECOMPRESSED}" = "1" ]; then
            rm -f \${read0} \${read1}
        fi
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================

workflow PRINSEQ_SUB {
    take:
    p
    reads

    main:
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = prinseq(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. エントリーポイント
// ==========================================

workflow PRINSEQ_ALL {
    p     = createNullParamsChannel()
    reads = createPairsChannel(params.prinseq_reads)

    out_ch = PRINSEQ_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

workflow {
    PRINSEQ_ALL()
}