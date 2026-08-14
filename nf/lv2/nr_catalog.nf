#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def RENAME_CATALOG_MAX_MEMORY = 4
def RENAME_CATALOG_MAX_THREADS = 1
def FILTER_FNA_MAX_MEMORY   = 8
def FILTER_FNA_MAX_THREADS  = 2

// 3. 上限値による動的クリッピング
params.nr_catalog_rename_catalog_memory  = Math.min(params.memory as Integer, RENAME_CATALOG_MAX_MEMORY)
params.nr_catalog_rename_catalog_threads = Math.min(params.threads as Integer, RENAME_CATALOG_MAX_THREADS)
params.nr_catalog_filter_fna_memory    = Math.min(params.memory as Integer, FILTER_FNA_MAX_MEMORY)
params.nr_catalog_filter_fna_threads   = Math.min(params.threads as Integer, FILTER_FNA_MAX_THREADS)

include { createNullParamsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { BUILD_REF_DB_PROT_SUB; CLUSTER_PROT_SUB } from "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// ==========================================
// 1. プロセス定義
// ==========================================

process rename_catalog {
    tag "rename_catalog"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.nr_catalog_rename_catalog_memory}"
    def threads = "${params.nr_catalog_rename_catalog_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(id), path(rep_faa)
    output:
        tuple val("nr_catalog"), path("nr_catalog.faa"), path("id_mapping.tsv")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        # FASTAのヘッダーを連番(NRCAT0000001...)に書き換えつつ、旧IDとの対応表(id_mapping.tsv)を作成するPythonスクリプト
        python3 - << 'EOF'
        import sys

        input_faa = "${rep_faa}"
        output_faa = "nr_catalog.faa"
        mapping_tsv = "id_mapping.tsv"

        count = 1
        with open(input_faa, "r") as fin, open(output_faa, "w") as fout, open(mapping_tsv, "w") as fmap:
            for line in fin:
                if line.startswith(">"):
                    old_id = line.strip()[1:].split()[0]
                    new_id = f"NRCAT{count:07d}"
                    fout.write(f">{new_id}\\n")
                    fmap.write(f"{new_id}\\t{old_id}\\n")
                    count += 1
                else:
                    fout.write(line)
        EOF
        """
}

process filter_fna_by_faa {
    tag "filter_fna"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.nr_catalog_filter_fna_memory}"
    def threads = "${params.nr_catalog_filter_fna_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(id), path(merged_fna)
        path(id_mapping) // id_mapping.tsv から旧IDを抽出して seqkit grep に使う
    output:
        tuple val("nr_catalog"), path("nr_catalog.fna")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        # mappingファイルから旧IDの列（2列目）を抽出して一時ファイルに落とし、seqkitに渡す
        awk '{print \$2}' ${id_mapping} > old_rep_ids.txt
        seqkit grep -f old_rep_ids.txt ${merged_fna} > nr_catalog_raw.fna
        
        # 塩基配列側もアミノ酸側と同様に連番ID（NRCAT...）へヘッダーを置換する
        python3 - << 'EOF'
        mapping = {}
        with open("${id_mapping}") as f:
            for line in f:
                new_id, old_id = line.strip().split("\\t")
                mapping[old_id] = new_id

        input_fna = "nr_catalog_raw.fna"
        output_fna = "nr_catalog.fna"

        with open(input_fna, "r") as fin, open(output_fna, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    orig_header = line.strip()[1:]
                    # seqkit grep はヒットした配列のヘッダーをそのまま返すため、先頭の単語（ID）でマッピングを引く
                    orig_id = orig_header.split()[0]
                    new_id = mapping.get(orig_id, orig_id)
                    fout.write(f">{new_id}\\n")
                else:
                    fout.write(line)
        EOF
        """
}

// ==========================================
// 2. サブワークフロー
// ==========================================

workflow NR_CATALOG_SUB {
    take:
    p
    merged_faa // [ "merged_all_samples", merged.faa ]
    merged_fna // [ "merged_all_samples", merged.fna ]

    main:
    // 1. タンパク質ベースのクラスタリング用 DB 作成 & クラスタリング実行
    cluster_in = merged_faa.map { id, path -> tuple("nr_prot", path) }
    ref_db     = BUILD_REF_DB_PROT_SUB(p, cluster_in)
    clustered  = CLUSTER_PROT_SUB(p, ref_db)

    // 2. 代表配列のIDリネーム & 対応表作成 (NRCAT0000001形式へ)
    // 出力: [ "nr_catalog", nr_catalog.faa, id_mapping.tsv ]
    renamed_catalog = rename_catalog(clustered)

    // 3. 塩基配列のフィルタリングとヘッダーの連番置換
    // renamed_catalog の中から [id, faa, mapping_tsv] を受け取り、mapping_tsv をフィルタに利用
    mapping_ch = renamed_catalog.map { id, faa, mapping -> tuple(id, mapping) }
    
    nr_fna = filter_fna_by_faa(merged_fna, mapping_ch.map { it[1] })

    emit:
    rep_faa = renamed_catalog.map { id, faa, mapping -> tuple(id, faa) }
    rep_fna = nr_fna
    mapping = mapping_ch.map { id, mapping -> tuple(id, mapping) }
}

// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow NR_CATALOG_ALL {
    p         = createNullParamsChannel()
    catalog_faa = createSeqsChannel(params.nr_catalog_faa)
    catalog_fna = createSeqsChannel(params.nr_catalog_fna)

    out_ch = NR_CATALOG_SUB(p, catalog_faa, catalog_fna)

    out_ch.rep_faa.view { i -> "REP_FAA: $i" }
    out_ch.rep_fna.view { i -> "REP_FNA: $i" }
    out_ch.mapping.view { i -> "MAPPING: $i" }
}

workflow {
    NR_CATALOG_ALL()
}