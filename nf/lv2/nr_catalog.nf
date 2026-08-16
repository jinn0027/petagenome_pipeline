#!/usr/bin/env nextflow
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

include { createNullParamsChannel; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { BUILD_REF_DB_PROT_SUB; CLUSTER_PROT_SUB } from "${params.petagenomeDir}/nf/lv1/mmseqs2.nf"

// ==========================================
// 1. プロセス定義
// ==========================================

process sanitize_and_rename {
    tag "sanitize_${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb     = "${params.nr_catalog_sanitize_memory}"
    def threads = "${params.nr_catalog_sanitize_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(id), path(faa), path(fna)
    output:
        tuple val(id), path("sanitized.faa"), path("sanitized.fna"), path("initial_mapping.tsv")
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

    def gb     = "${params.nr_catalog_filter_fna_memory}"
    def threads = "${params.nr_catalog_filter_fna_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        path sanitized_fnas       // 複数のfnaファイルをまとめて受け取る
        tuple val(rep_id), path(rep_faa)
    output:
        tuple val("nr_catalog"), path("nr_catalog.fna")
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
rep_ids = set()
with open("${rep_faa}", "r") as f:
    for line in f:
        if line.startswith(">"):
            rep_id = line.strip()[1:].split()[0]
            rep_ids.add(rep_id)

output_fna = "nr_catalog.fna"

with open(output_fna, "w") as fout:
    for input_fna in "${sanitized_fnas}".split():
        with open(input_fna, "r") as fin:
            write_flag = False
            for line in fin:
                if line.startswith(">"):
                    orig_id = line.strip()[1:].split()[0]
                    if orig_id in rep_ids:
                        fout.write(line)
                        write_flag = True
                    else:
                        write_flag = False
                else:
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
    input_ch = faa.join(fna)
    sanitized = sanitize_and_rename(input_ch) // [id, sanitized_faa, sanitized_fna, mapping]

    // 1. MMseqs2 用に入力を成形
    cluster_in = sanitized.map { id, faa_p, fna_p, mapping -> tuple("nr_prot", faa_p) }
    ref_db     = BUILD_REF_DB_PROT_SUB(p, cluster_in)
    clustered_raw = CLUSTER_PROT_SUB(p, ref_db) 

    // どのような形式で返ってきても確実に [id, rep_faa] の2要素タプルに整形する
    clustered = clustered_raw.map { item ->
        if (item instanceof List) {
            def path_val = item.find { it instanceof java.nio.file.Path || it instanceof File || (it.toString().endsWith('.fasta') || it.toString().endsWith('.faa')) }
            def id_val = item[0] instanceof CharSequence ? item[0] : "nr_prot"
            return tuple(id_val, path_val)
        }
        return item
    }

    // 2. 全サンプルの sanitized_fna をリストとして収集して渡す
    sanitized_fna_list = sanitized.map { id, faa_p, fna_p, mapping -> fna_p }.collect()

    nr_fna = filter_fna_by_faa(sanitized_fna_list, clustered)

    mapping_ch = sanitized.map { id, faa_p, fna_p, mapping -> tuple(id, mapping) }

    emit:
    rep_faa = clustered
    rep_fna = nr_fna
    mapping = mapping_ch
}

// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow NR_CATALOG_ALL {
    p             = createNullParamsChannel()
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