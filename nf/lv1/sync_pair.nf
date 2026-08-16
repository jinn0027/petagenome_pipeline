#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. このモジュール・タスク固有の推奨・上限値
def SYNC_MAX_MEMORY  = 32  // 辞書やハッシュを持つため少し余裕を持たせる
def SYNC_MAX_THREADS = 8   

// 3. 上限値による動的クリッピング
params.sync_pair_memory  = Math.min(params.memory as Integer, SYNC_MAX_MEMORY)
params.sync_pair_threads = Math.min(params.threads as Integer, SYNC_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process sync_pair {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    
    def gb = "${params.sync_pair_memory}"
    def threads = "${params.sync_pair_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("${pair_id}/sync_{1,2}.fastq.gz", arity: '2')

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # Pythonによる厳密なペア同期処理
        python3 - << 'EOF'
        import gzip

        def clean_qname(title):
            # ヘッダーの最初の単語（スペース前）を取り出し、末尾の /1 または /2 のみを除去
            first_word = title.strip().split()[0]
            if first_word.endswith('/1') or first_word.endswith('/2'):
                return first_word[:-2]
            return first_word

        r1_path = "${reads[0]}"
        r2_path = "${reads[1]}"
        out1_path = "${pair_id}/sync_1.fastq.gz"
        out2_path = "${pair_id}/sync_2.fastq.gz"

        print("Indexing R2 reads...")
        r2_dict = {}
        with gzip.open(r2_path, 'rt') as f2:
            while True:
                l1 = f2.readline()
                if not l1: break
                l2 = f2.readline()
                l3 = f2.readline()
                l4 = f2.readline()
                
                qname = clean_qname(l1)
                r2_dict[qname] = (l1, l2, l3, l4)

        print(f"Total R2 records indexed: {len(r2_dict)}")

        print("Filtering and writing synced pairs...")
        matched_count = 0
        with gzip.open(r1_path, 'rt') as f1, \
             gzip.open(out1_path, 'wt') as out1, \
             gzip.open(out2_path, 'wt') as out2:
            
            while True:
                l1 = f1.readline()
                if not l1: break
                l2 = f1.readline()
                l3 = f1.readline()
                l4 = f1.readline()
                
                qname = clean_qname(l1)
                
                # R1とR2の両方に存在する（ペアが揃っている）場合のみ出力
                if qname in r2_dict:
                    r2_rec = r2_dict[qname]
                    
                    # R1を書く
                    out1.write(l1)
                    out1.write(l2)
                    out1.write(l3)
                    out1.write(l4)
                    
                    # 対応するR2を書く
                    out2.write(r2_rec[0])
                    out2.write(r2_rec[1])
                    out2.write(r2_rec[2])
                    out2.write(r2_rec[3])
                    
                    matched_count += 1

        print(f"Successfully synced {matched_count} pairs.")
        EOF
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow SYNC_PAIR_SUB {
    take:
    p
    reads

    main:
    in_ch = p.combine(reads).map { p_val, pair_id, reads_path ->
        tuple(p_val, pair_id, reads_path)
    }

    out = sync_pair(in_ch)

    emit:
    out = out
}

// ==========================================
// 2. エントリーポイント (-entry)
// ==========================================
workflow SYNC_PAIR_ALL {
    p      = createNullParamsChannel()
    reads = createPairsChannel(params.sync_pair_reads)

    out_ch = SYNC_PAIR_SUB(p, reads)
    out_ch.out.view { i -> "$i" }
}

workflow {
    SYNC_PAIR_ALL()
}