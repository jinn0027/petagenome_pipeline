#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_GG2_WGS_MAX_MEMORY = 64
def QIIME2_GG2_WGS_MAX_THREADS = 32

params.qiime2_gg2_wgs_memory  = Math.min(params.memory as Integer, QIIME2_GG2_WGS_MAX_MEMORY)
params.qiime2_gg2_wgs_threads = Math.min(params.threads as Integer, QIIME2_GG2_WGS_MAX_THREADS)

// ミスマッチ許容数のデフォルト設定（-1: 制限なし, 0: 完全一致, 1: 1塩基違いまで許可 ...）
if (!params.containsKey('qiime2_gg2_wgs_max_mismatch')) {
    params.qiime2_gg2_wgs_max_mismatch = -1
}

// リファレンスファイルのパスのデフォルト設定
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy     = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process qiime2_greengenes2_wgs {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif" // python, pandas等を利用
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_gg2_wgs_memory}"
    def threads = "${params.qiime2_gg2_wgs_threads}"
    def max_mismatch = params.qiime2_gg2_wgs_max_mismatch
    
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(alignment_file) // Bowtie2などのアライメント結果(SAM/BAM)
        path backbone_fna
        path taxonomy

    output:
        tuple val(pair_id), 
              path("${pair_id}/feature-table.tsv"), 
              path("${pair_id}/taxonomy.tsv"),
              path("${pair_id}/taxonomy_counts.tsv") // ← 追加

    script:
        """
        # Pythonの非推奨警告を抑制
        export PYTHONWARNINGS="ignore"

        # コンテナ特有のキャッシュ・ホームディレクトリ競合を防ぐ環境変数
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        # 1. .qza 形式のタクソノミーファイルから TSV をエクスポート（必要に応じて）
        qiime tools export \
            --input-path ${taxonomy} \
            --output-path exported_taxonomy

        TAX_TSV=""
        if [ -f exported_taxonomy/taxonomy.tsv ]; then
            TAX_TSV="exported_taxonomy/taxonomy.tsv"
        elif [ -f exported_taxonomy/consensus_assignments.tsv ]; then
            TAX_TSV="exported_taxonomy/consensus_assignments.tsv"
        else
            TAX_TSV=\$(find exported_taxonomy -name "*.tsv" | head -n 1)
        fi

        # 2. Pythonスクリプトにより、アライメント結果をパースし、
        #    ミスマッチ数制限（max_mismatch）に応じてフィルタリングして集計する
        python3 - <<EOF
        import pandas as pd

        max_mismatch = ${max_mismatch}
        print(f"Loading taxonomy from \${TAX_TSV}...")
        print(f"Mismatch threshold: {max_mismatch} (-1 means no restriction)")

        tax_df = pd.read_csv("\${TAX_TSV}", sep="\\t", header=None, index_col=0)
        tax_dict = tax_df[1].to_dict()

        counts = {}
        mapped = 0
        unmapped = 0

        print("Parsing alignment file (${alignment_file})...")
        with open("${alignment_file}", "r") as f:
            for line in f:
                if line.startswith("@"):
                    continue
                parts = line.strip().split("\\t")
                if len(parts) < 3:
                    continue
                
                flag = int(parts[1])
                ref_id = parts[2] # ヒットしたバックボーン配列のID

                # 未マッピング (Bit 4: 0x4) または "*" の場合はスキップ
                if (flag & 4) or ref_id == "*":
                    unmapped += 1
                    continue

                # ミスマッチ数 (NM:i:N) の判定
                if max_mismatch >= 0:
                    is_valid_match = False
                    for tag in parts[11:]:
                        if tag.startswith("NM:i:"):
                            try:
                                nm_val = int(tag.split(":")[2])
                                if nm_val <= max_mismatch:
                                    is_valid_match = True
                            except ValueError:
                                pass
                            break
                    
                    if not is_valid_match:
                        unmapped += 1
                        continue

                mapped += 1
                taxon = tax_dict.get(ref_id, "k__Unassigned; p__; c__; o__; f__; g__; s__")
                counts[taxon] = counts.get(taxon, 0) + 1

        print(f"Accepted reads: {mapped}, Filtered/Unmapped reads: {unmapped}")

        taxa_list = list(counts.keys())
        freq_list = list(counts.values())

        # taxonomy.tsv の出力
        tax_out = pd.DataFrame({
            "Feature ID": taxa_list,
            "Taxon": taxa_list
        })
        tax_out.to_csv("${pair_id}/taxonomy.tsv", sep="\\t", index=False)

        # feature-table.tsv の出力
        table_out = pd.DataFrame({
            "Taxonomy": taxa_list,
            "${pair_id}": freq_list
        })
        table_out.to_csv("${pair_id}/feature-table.tsv", sep="\\t", index=False)

        # 【追加】カウント数の大きい順（降順）にソートした taxonomy_counts.tsv の出力
        summary_out = pd.DataFrame({
            "Taxonomy": taxa_list,
            "${pair_id}": freq_list
        }).sort_values(by="${pair_id}", ascending=False)
        
        summary_out.to_csv("${pair_id}/taxonomy_counts.tsv", sep="\\t", index=False)

        print("Done successfully.")
        EOF
        """
}

// ==========================================
// 1. サブワークフロー
// ==========================================
workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    alignment_ch
    backbone_fna
    taxonomy

    main:
    out = qiime2_greengenes2_wgs(
        alignment_ch,
        backbone_fna,
        taxonomy
    )

    emit:
    out = out
}

// ==========================================
// 2. コマンドライン用エントリーポイント
// ==========================================
workflow QIIME2_GREENGENES2_WGS_ALL {
    p = createNullParamsChannel()
    
    backbone_ch = Channel.value(file(params.qiime2_gg2_wgs_backbone_fna, checkIfExists: true))
    taxonomy_ch = Channel.value(file(params.qiime2_gg2_wgs_taxonomy, checkIfExists: true))

    alignment_dummy_ch = Channel.fromPath("${params.output}/MAP_SUB/*", checkIfExists: false)
        .map { dir_path ->
            def pair_id = dir_path.name
            def sam_file = file("${dir_path}/*.sam", checkIfExists: false)
            return tuple(pair_id, sam_file)
        }

    QIIME2_GREENGENES2_WGS_SUB(
        p,
        alignment_dummy_ch,
        backbone_ch,
        taxonomy_ch
    )
}

workflow {
    QIIME2_GREENGENES2_WGS_ALL()
}