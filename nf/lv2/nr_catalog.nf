#!/usr/init/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def SANITIZE_MAX_MEMORY      = 8
def SANITIZE_MAX_THREADS     = 2
def FILTER_FNA_MAX_MEMORY    = 8
def FILTER_FNA_MAX_THREADS   = 2

// 3. 上限値による動的クリッピング
params.nr_catalog_sanitize_memory    = Math.min(params.memory as Integer, SANITIZE_MAX_MEMORY)
params.nr_catalog_sanitize_threads   = Math.min(params.threads as Integer, SANITIZE_MAX_THREADS)
params.nr_catalog_filter_fna_memory  = Math.min(params.memory as Integer, FILTER_FNA_MAX_MEMORY)
params.nr_catalog_filter_fna_threads = Math.min(params.threads as Integer, FILTER_FNA_MAX_THREADS)

include { createNullParamsChannel; createSeqsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { BUILD_REF_DB_PROT_SUB; CLUSTER_PROT_SUB } from "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// ==========================================
// 1. プロセス定義
// ==========================================

// 複数ファイルの単純な結合用プロセス
process merge_input_sequences {
    tag "merge_inputs"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    def gb        = "${params.nr_catalog_sanitize_memory}"
    def threads = "${params.nr_catalog_sanitize_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads

    input:
        path(faas)
        path(fnas)
    output:
        tuple path("merged.faa"), path("merged.fna")
    script:
        """
        cat ${faas} > merged.faa
        cat ${fnas} > merged.fna
        """
}

process sanitize_and_rename {
    tag "sanitize_all"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb        = "${params.nr_catalog_sanitize_memory}"
    def threads = "${params.nr_catalog_sanitize_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple path(faa), path(fna)
    output:
        tuple path("sanitized.faa"), path("sanitized.fna"), path("initial_mapping.tsv")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
input_faa = "${faa}"
input_fna = "${fna}"
out_faa = "sanitized.faa"
out_fna = "sanitized.fna"
mapping = "initial_mapping.tsv"

faa_dict = {}
with open(input_faa, "r") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header:
                faa_dict[header.strip()[1:].split()[0]] = "".join(seq_lines)
            header = line
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header:
        faa_dict[header.strip()[1:].split()[0]] = "".join(seq_lines)

fna_dict = {}
with open(input_fna, "r") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header:
                fna_dict[header.strip()[1:].split()[0]] = "".join(seq_lines)
            header = line
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header:
        fna_dict[header.strip()[1:].split()[0]] = "".join(seq_lines)

common_keys = [k for k in faa_dict.keys() if k in fna_dict]

count = 1
with open(out_faa, "w") as ofaa, open(out_fna, "w") as ofna, open(mapping, "w") as fmap:
    for old_id in common_keys:
        new_id = f"GENE{count:09d}"
        ofaa.write(f">{new_id}\\n{faa_dict[old_id]}\\n")
        ofna.write(f">{new_id}\\n{fna_dict[old_id]}\\n")
        fmap.write(f"{new_id}\\t{old_id}\\n")
        count += 1
EOF
        """
}

process filter_fna_by_faa {
    tag "filter_fna"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb        = "${params.nr_catalog_filter_fna_memory}"
    def threads = "${params.nr_catalog_filter_fna_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple path(sanitized_fna), val(rep_id), path(rep_faa)

    output:
        tuple val("nr_catalog"), path("nr_catalog.fna")

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
valid_ids = set()
with open("${rep_faa}", "r") as f:
    for line in f:
        if line.startswith(">"):
            valid_ids.add(line.strip()[1:].split()[0])

with open("${sanitized_fna}", "r") as fin, open("nr_catalog.fna", "w") as fout:
    write_flag = False
    for line in fin:
        if line.startswith(">"):
            seq_id = line.strip()[1:].split()[0]
            write_flag = (seq_id in valid_ids)
        if write_flag:
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
    faa
    fna

    main:
    // 1. 全てのFAA/FNAを収集してまずマージ
    all_faa = faa.map { id, path -> path }.collect()
    all_fna = fna.map { id, path -> path }.collect()
    
    merged = merge_input_sequences(all_faa, all_fna)
    merged.view { "DEBUG [NR_CATALOG merged]: $it" }

    // 2. マージされたファイルに対して1回だけサニタイズを実行
    sanitized = sanitize_and_rename(merged)
    sanitized.view { "DEBUG [NR_CATALOG sanitized]: $it" }

    // 複数の処理（cluster, filter, mapping）で使うため、チャンネルを事前に分岐
    sanitized_for_cluster = sanitized.map { faa, fna, map -> tuple(faa, fna, map) }
    sanitized_for_filter  = sanitized.map { faa, fna, map -> tuple(faa, fna, map) }
    sanitized_for_mapping = sanitized.map { faa, fna, map -> tuple(faa, fna, map) }

    // 3. MMseqs2 用の入力: [ "nr_prot", sanitized.faa ]
    cluster_in = sanitized_for_cluster.map { faa_p, fna_p, mapping -> tuple("nr_prot", faa_p) }
    cluster_in.view { "DEBUG [NR_CATALOG cluster_in]: $it" }

    // 4. クラスタリング用のDB作成
    ref_db = BUILD_REF_DB_PROT_SUB(p, cluster_in)
    clustered_faa = CLUSTER_PROT_SUB(p, ref_db).out.map { id, fasta, tsv -> tuple(id, fasta) }
    clustered_faa.view { "DEBUG [NR_CATALOG clustered_faa]: $it" }

    // 5. filter_in の構築および filter_fna_by_faa 呼び出し
    // .combine() は要素をフラット化するため、5つの要素を直接受け取る
    filter_in = sanitized_for_filter
        .combine(clustered_faa)
        .map { sanitized_faa, sanitized_fna, mapping, rep_id, rep_faa ->
            tuple(sanitized_fna, rep_id, rep_faa)
        }
    filter_in.view { "DEBUG [NR_CATALOG filter_in]: $it" }

    nr_fna = filter_fna_by_faa(filter_in)
    nr_fna.view { "DEBUG [NR_CATALOG nr_fna]: $it" }

    // 6. マッピングの出力
    nr_mapping = sanitized_for_mapping.map { faa_p, fna_p, mapping -> tuple("nr_catalog", mapping) }
    nr_mapping.view { "DEBUG [NR_CATALOG nr_mapping]: $it" }
    
    emit:
    rep_faa = clustered_faa      // [ "nr_prot", Path ]
    rep_fna = nr_fna               // [ "nr_catalog", Path ]
    mapping = nr_mapping           // [ "nr_catalog", Path ]
}

workflow NR_CATALOG_ALL {
    p                = createNullParamsChannel()
    catalog_faa = createSeqsChannel(params.nr_catalog_faa)
    catalog_fna = createSeqsChannel(params.nr_catalog_fna)

    out_ch = NR_CATALOG_SUB(p, catalog_faa, catalog_fna)

    out_ch.rep_faa.view { "DEBUG [REP_FAA]: $it" }
    out_ch.rep_fna.view { "DEBUG [REP_FNA]: $it" }
    out_ch.mapping.view { "DEBUG [MAPPING]: $it" }
}

workflow {
    NR_CATALOG_ALL()
}