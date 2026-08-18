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
    def gb      = "${params.nr_catalog_sanitize_memory}"
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

    def gb      = "${params.nr_catalog_sanitize_memory}"
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

    def gb      = "${params.nr_catalog_filter_fna_memory}"
    def threads = "${params.nr_catalog_filter_fna_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple path(sanitized_fnas), val(rep_id), path(rep_faa)

    output:
        tuple val("nr_catalog"), path("nr_catalog.fna")

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
import glob, os

rep_ids = set()
with open("${rep_faa}", "r") as f:
    for line in f:
        if line.startswith(">"):
            rid = line.strip()[1:].split()[0]
            rep_ids.add(rid)

output_fna = "nr_catalog.fna"

input_fna_files = [f for f in glob.glob("*.fna") if f != output_fna]

with open(output_fna, "w") as fout:
    for input_fna in sorted(input_fna_files):
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
    // 1. 全てのFAA/FNAを収集してまずマージ
    all_faa = faa.map { id, path -> path }.collect()
    all_fna = fna.map { id, path -> path }.collect()
    
    merged = merge_input_sequences(all_faa, all_fna)

    // 2. マージされたファイルに対して1回だけサニタイズを実行
    sanitized = sanitize_and_rename(merged)

    // 3. MMseqs2 用の入力: [ "nr_prot", sanitized.faa ]
    cluster_in = sanitized.map { faa_p, fna_p, mapping -> tuple("nr_prot", faa_p) }

    // 4. クラスタリング用のDB作成
    ref_db = BUILD_REF_DB_PROT_SUB(p, cluster_in)
    clustered_faa = CLUSTER_PROT_SUB(p, ref_db).out.map { id, fasta, tsv -> tuple(id, fasta) }

    // 4. フィルタリング用の入力組み立て ([sanitized.fna], rep_id, rep_faa)
    filter_in = sanitized.map { faa_p, fna_p, mapping -> [fna_p] }.combine(clustered_faa)
    
    nr_fna = filter_fna_by_faa(filter_in)

    // 5. マッピングの出力（1つになった initial_mapping.tsv をそのまま利用）
    nr_mapping = sanitized.map { faa_p, fna_p, mapping -> tuple("nr_catalog", mapping) }
    
    emit:
    rep_faa = clustered_faa      // [ "nr_prot", Path ]
    rep_fna = nr_fna             // [ "nr_catalog", Path ]
    mapping = nr_mapping         // [ "nr_catalog", Path ]
}

workflow NR_CATALOG_ALL {
    p                = createNullParamsChannel()
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