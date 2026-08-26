#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================
// 0. グローバル関数定義
// ==========================================
def checkMandatoryFile(paramName, pathStr) {
    if (!pathStr) { error "【パラメータ不足】 '--${paramName}' が指定されていません。" }
    def f = file(pathStr)
    if (!f.exists()) { error "【ファイル非存在エラー】 '--${paramName}' で指定されたファイルが存在しません: ${pathStr}" }
    return f
}

// 1. 全体デフォルト値の定義
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def SYNC_PAIR_MAX_MEMORY = 32
def SYNC_PAIR_MAX_THREADS = 8
def MERGE_FASTA_MAX_MEMORY  = 16
def MERGE_FASTA_MAX_THREADS = 4
def FILTER_CATALOG_MAX_MEMORY  = 8
def FILTER_CATALOG_MAX_THREADS = 2

// 3. 上限値による動的クリッピング
params.sync_pair_memory     = Math.min(params.memory as Integer, SYNC_PAIR_MAX_MEMORY)
params.sync_pair_threads    = Math.min(params.threads as Integer, SYNC_PAIR_MAX_THREADS)
params.merge_fasta_memory   = Math.min(params.memory as Integer, MERGE_FASTA_MAX_MEMORY)
params.merge_fasta_threads  = Math.min(params.threads as Integer, MERGE_FASTA_MAX_THREADS)
params.filter_catalog_memory    = Math.min(params.memory as Integer, FILTER_CATALOG_MAX_MEMORY)
params.filter_catalog_threads   = Math.min(params.threads as Integer, FILTER_CATALOG_MAX_THREADS)

include { createNullParamsChannel; createSeqsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// 各種サブワークフローのインポート
include { FASTP_SUB }                  from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB }          from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { ASSEMBLY_SUB }             from "${params.petagenomeDir}/nf/lv2/assembly.nf"
include { PRODIGAL_SUB }             from "${params.petagenomeDir}/nf/lv1/prodigal.nf"
include { ANNOTATE_CATALOG_SUB }     from "${params.petagenomeDir}/nf/lv2/annotate_catalog_prot.nf"
include { NR_CATALOG_SUB }             from "${params.petagenomeDir}/nf/lv2/nr_catalog.nf"
include { ALIGN_SAMPLES_TO_CATALOG_SUB } from "${params.petagenomeDir}/nf/lv2/align_samples_to_catalog.nf"
include { ANNOTATE_SAMPLES_SUB }       from "${params.petagenomeDir}/nf/lv2/annotate_samples.nf"
include { SUMMARIZE_KO_TAXID_SUB }     from "${params.petagenomeDir}/nf/lv2/summarize_ko_taxid.nf"

// プロセス定義
process merge_fasta {
    tag "${ext}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.merge_fasta_memory}"; def threads = "${params.merge_fasta_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"; cpus params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input: tuple path(fasta_files, stageAs: "?/*"), val(ext)
    output: tuple val("merged_all_samples"), path("combined_orfs.${ext}")
    script: """ echo "${processProfile(task)}" | tee prof.txt; cat ${fasta_files} > combined_orfs.${ext} """
}

// ペアエンドリードを1つに結合するプロセス（MMseqs2のクエリ用）
process cat_paired_reads {
    tag "${qry_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    input: tuple val(qry_id), path(reads) // reads は [r1, r2] のリスト
    output: tuple val(qry_id), path("combined.fastq.gz")
    script:
    """
    cat ${reads[0]} ${reads[1]} > combined.fastq.gz
    """
}

process filter_faa_catalog_by_hits {
    tag "faa"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.filter_catalog_memory}"; def threads = "${params.filter_catalog_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"; cpus params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input: tuple path(fasta_file), val(ext), path(hit_ids)
    output: tuple val("filtered_catalog"), path("filtered_catalog_faa.faa")
    script: """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        valid_ids = {line.strip() for line in open("${hit_ids}") if line.strip()}
        with open("${fasta_file}", "r") as fin, open("filtered_catalog_faa.faa", "w") as fout:
            write_flag = False
            for line in fin:
                if line.startswith(">"): write_flag = (line.strip()[1:].split()[0] in valid_ids)
                if write_flag: fout.write(line)
        EOF
    """
}

process filter_fna_catalog_by_hits {
    tag "fna"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.filter_catalog_memory}"; def threads = "${params.filter_catalog_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"; cpus params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input: tuple path(fasta_file), val(ext), path(hit_ids)
    output: tuple val("filtered_catalog"), path("filtered_catalog_fna.fna")
    script: """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        valid_ids = {line.strip() for line in open("${hit_ids}") if line.strip()}
        with open("${fasta_file}", "r") as fin, open("filtered_catalog_fna.fna", "w") as fout:
            write_flag = False
            for line in fin:
                if line.startswith(">"): write_flag = (line.strip()[1:].split()[0] in valid_ids)
                if write_flag: fout.write(line)
        EOF
    """
}

process filter_mapping_by_hits {
    tag "filter_mapping"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.filter_catalog_memory}"; def threads = "${params.filter_catalog_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"; cpus params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input: tuple path(id_mapping), path(hit_ids)
    output: tuple val("filtered_mapping"), path("filtered_id_mapping.tsv")
    script: """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        valid_ids = {line.strip() for line in open("${hit_ids}") if line.strip()}
        with open("${id_mapping}", "r") as fin, open("filtered_id_mapping.tsv", "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if parts and parts[0] in valid_ids: fout.write(line)
        EOF
    """
}

// ==========================================
// サブワークフロー定義
// ==========================================

// (0-1) 各サンプルの前処理・アセンブリ・ORF予測（ORFS解析用）
workflow PREPARE_ORFS_SUB {
    take: p; host_ref; reads

    main:
    fp = FASTP_SUB(p, reads)
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)
    asm_res = ASSEMBLY_SUB(p, host_removed.reads)
    orf_res = PRODIGAL_SUB(p, asm_res.flt_seqs)

    emit:
    raw_reads          = reads
    fastp_reads        = fp.out
    host_removed_reads = host_removed.reads
    contigs            = asm_res.asm
    flt_seqs           = asm_res.flt_seqs
    orfs               = orf_res.out
}

// (0-2) クリーニング済みリードの準備（READS解析用）
workflow PREPARE_CLEAN_READS_SUB {
    take: p; host_ref; reads

    main:
    fp = FASTP_SUB(p, reads)
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)

    emit:
    reads = host_removed.reads
}

// (1) カタログ構築ワークフロー
workflow BUILD_CATALOG_SUB {
    take: p; ref_or_db; taxid_map; ko_map; orf_res_out

    main:
    sample_orf_fna_ch = orf_res_out.map { qry_id, faa, fna, gbk -> tuple(qry_id, fna) }
    all_faa_files = orf_res_out.map { qry_id, faa, fna, gbk -> faa }.collect()
    all_fna_files = orf_res_out.map { qry_id, faa, fna, gbk -> fna }.collect()

    merge_inputs_ch = all_faa_files.map { tuple(it, "faa") }.mix(all_fna_files.map { tuple(it, "fna") })
    merged_fastas = merge_fasta(merge_inputs_ch)

    merged_faa_fasta = merged_fastas.filter { id, path -> path.name.endsWith('.faa') }
    merged_fna_fasta = merged_fastas.filter { id, path -> path.name.endsWith('.fna') }

    nr_catalog_res = NR_CATALOG_SUB(p, merged_faa_fasta, merged_fna_fasta)
    nr_rep_faa_path = nr_catalog_res.rep_faa.map { id, path -> path }
    nr_rep_fna_path = nr_catalog_res.rep_fna.map { id, path -> path }

    annotation_res = ANNOTATE_CATALOG_SUB(p, ref_or_db, nr_catalog_res.rep_faa, taxid_map, ko_map)
    annotated_path_ch = annotation_res.annotated.map { id, path -> path }
    hit_ids_path_ch   = annotation_res.hit_ids.map { id, path -> path }

    faa_to_filter = nr_rep_faa_path.map { tuple(it, "faa") }.combine(hit_ids_path_ch)
    fna_to_filter = nr_rep_fna_path.map { tuple(it, "fna") }.combine(hit_ids_path_ch)

    filtered_rep_faa = filter_faa_catalog_by_hits(faa_to_filter)
    filtered_rep_fna = filter_fna_catalog_by_hits(fna_to_filter)

    mapping_path_ch = nr_catalog_res.mapping.map { id, path -> path }
    mapping_to_filter = mapping_path_ch.combine(hit_ids_path_ch)
    filtered_mapping_res = filter_mapping_by_hits(mapping_to_filter)

    emit:
    nr_catalog_faa     = nr_catalog_res.rep_faa
    nr_catalog_fna     = nr_catalog_res.rep_fna
    annotated          = annotated_path_ch
    hit_ids            = hit_ids_path_ch
    filtered_rep_faa   = filtered_rep_faa
    filtered_rep_fna   = filtered_rep_fna
    filtered_mapping   = filtered_mapping_res
    sample_orf_fna     = sample_orf_fna_ch
}

// (2) カタログへのマッピング・アノテーション・集計ワークフロー
workflow ANNOTATE_AND_SUMMARIZE_SAMPLES_SUB {
    take: p; filtered_rep_fna; target_reads_ch; annotated_path_ch; ko_name_map; taxid_name_map

    main:
    samples_mapping_res = ALIGN_SAMPLES_TO_CATALOG_SUB(p, filtered_rep_fna, target_reads_ch)
    
    // 【期待される形式】: tuple(val(qry_id), path(m8_file))
    samples_mapping_res.out.view { 
        "[DEBUG VIEW] samples_mapping_res.out\n" +
        "  - 期待される形式: tuple(val(qry_id), path(m8_file)) [例: ['01_sampleA', /path/to/out.m8]]\n" +
        "  - 実際の値      : ${it}\n" +
        "  - クラス        : ${it?.getClass()}"
    }

    // 【期待される形式】: path(annotated_tsv_file) または tuple(val(id), path(tsv)) 等
    annotated_path_ch.view { 
        "[DEBUG VIEW] annotated_path_ch\n" +
        "  - 期待される形式: path(annotated_tsv_file) [例: /path/to/nr_prot_annotated.tsv]\n" +
        "  - 実際の値      : ${it}\n" +
        "  - クラス        : ${it?.getClass()}"
    }
    
    samples_anno_res = ANNOTATE_SAMPLES_SUB(p, samples_mapping_res.out, annotated_path_ch)
    
    // 【期待される形式】: tuple(val(qry_id), path(annotated_sample_tsv))
    samples_anno_res.samples_anno.view {
        "[DEBUG VIEW] samples_anno_res.samples_anno\n" +
        "  - 期待される形式: tuple(val(qry_id), path(annotated_sample_tsv))\n" +
        "  - 実際の値      : ${it}\n" +
        "  - クラス        : ${it?.getClass()}"
    }

    summary_res = SUMMARIZE_KO_TAXID_SUB(p, samples_anno_res.samples_anno, ko_name_map, taxid_name_map)

    emit:
    samples_mapping    = samples_mapping_res.out
    samples_annotation = samples_anno_res.samples_anno
    ko_summary         = summary_res.ko_summary
    tax_summary        = summary_res.tax_summary
}


// ==========================================
// 実行用エントリポイントワークフロー
// ==========================================

def setupReadsAndRefs() {
    if (!params.bacteriome_pipeline_reads) { error "【パラメータ不足】 '--bacteriome_pipeline_reads' が指定されていません。" }
    if (!params.remove_host_ref_fasta_or_db) { error "【パラメータ不足】 '--remove_host_ref_fasta_or_db' が指定されていません。" }
    if (!params.annotate_catalog_prot_ref_or_db) { error "【パラメータ不足】 '--annotate_catalog_prot_ref_or_db' が指定されていません。" }

    def taxid_map      = checkMandatoryFile('taxid_map_path', params.taxid_map_path)
    def ko_map         = checkMandatoryFile('ko_map_path', params.ko_map_path)
    def ko_name_map    = checkMandatoryFile('ko_name_map_path', params.ko_name_map_path)
    def taxid_name_map = checkMandatoryFile('taxid_name_map_path', params.taxid_name_map_path)

    def p = createNullParamsChannel()
    def reads_list = params.bacteriome_pipeline_reads.split(';')
    def individual_channels = reads_list.collect { Channel.fromFilePairs(it, checkIfExists: true) }
    def reads_mixed = individual_channels.first()
    if (individual_channels.size() > 1) { 
        individual_channels.tail().each { reads_mixed = reads_mixed.mix(it) } 
    }
    
    def index = 0
    def reads = reads_mixed.map { id, pair -> tuple("${String.format('%02d', index++)}_${id}", pair) }

    def host_ref  = createSeqsChannel(params.remove_host_ref_fasta_or_db)
    def ref_or_db = createSeqsChannel(params.annotate_catalog_prot_ref_or_db)

    return [p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads]
}

// 統合実行 (1) + (2) : 連続実行
workflow BACTERIOME_PIPELINE_ALL {
    def (p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads) = setupReadsAndRefs()

    orfs_res = PREPARE_ORFS_SUB(p, host_ref, reads)
    catalog_res = BUILD_CATALOG_SUB(p, ref_or_db, taxid_map, ko_map, orfs_res.orfs)
    summary_res = ANNOTATE_AND_SUMMARIZE_SAMPLES_SUB(p, catalog_res.filtered_rep_fna, catalog_res.sample_orf_fna, catalog_res.annotated, ko_name_map, taxid_name_map)
}

// カタログ作成のみ (1)
workflow BACTERIOME_PIPELINE_BUILD_CATALOG {
    def (p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads) = setupReadsAndRefs()

    orfs_res = PREPARE_ORFS_SUB(p, host_ref, reads)
    BUILD_CATALOG_SUB(p, ref_or_db, taxid_map, ko_map, orfs_res.orfs)
}

// サンプル解析（ORF配列ベース）
workflow BACTERIOME_PIPELINE_ANALYZE_ORFS {
    def (p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads) = setupReadsAndRefs()

    def catalog_fna_path        = checkMandatoryFile('catalog_fna', params.catalog_fna)
    def catalog_annotated_path = checkMandatoryFile('catalog_annotated', params.catalog_annotated)

    filtered_rep_fna_ch = Channel.value(tuple("filtered_catalog", file(catalog_fna_path)))
    annotated_ch         = Channel.value(file(catalog_annotated_path))

    orfs_res = PREPARE_ORFS_SUB(p, host_ref, reads)
    sample_orf_fna_ch = orfs_res.orfs.map { qry_id, faa, fna, gbk -> tuple(qry_id, fna) }

    ANNOTATE_AND_SUMMARIZE_SAMPLES_SUB(p, filtered_rep_fna_ch, sample_orf_fna_ch, annotated_ch, ko_name_map, taxid_name_map)
}

// サンプル解析（クリーニング済みリード直接ベース・アセンブリ/ORF予測なし）
workflow BACTERIOME_PIPELINE_ANALYZE_READS {
    def (p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads) = setupReadsAndRefs()

    def catalog_fna_path        = checkMandatoryFile('catalog_fna', params.catalog_fna)
    def catalog_annotated_path = checkMandatoryFile('catalog_annotated', params.catalog_annotated)

    filtered_rep_fna_ch = Channel.value(tuple("filtered_catalog", file(catalog_fna_path)))
    annotated_ch         = Channel.value(file(catalog_annotated_path))

    // 前処理（fastp + ホスト除去）のみを実行
    clean_reads_res = PREPARE_CLEAN_READS_SUB(p, host_ref, reads)

    // ペアエンドのリード（R1/R2）を1つのファイルに結合
    combined_reads_ch = cat_paired_reads(clean_reads_res.reads)

    // 結合されたクリーンリードを直接マッピングモジュールへ投入
    ANNOTATE_AND_SUMMARIZE_SAMPLES_SUB(p, filtered_rep_fna_ch, combined_reads_ch, annotated_ch, ko_name_map, taxid_name_map)
}

workflow { BACTERIOME_PIPELINE_ALL() }