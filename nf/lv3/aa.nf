#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
include { ASSEMBLY_SUB }               from "${params.petagenomeDir}/nf/lv2/assembly.nf"
include { PRODIGAL_SUB }               from "${params.petagenomeDir}/nf/lv1/prodigal.nf"
include { ANNOTATE_CATALOG_SUB }       from "${params.petagenomeDir}/nf/lv2/annotate_catalog_prot.nf"
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

process filter_catalog_by_hits {
    tag "${ext}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output
    def gb = "${params.filter_catalog_memory}"; def threads = "${params.filter_catalog_threads}"
    memory params.executor == "sge" ? null : "${gb} GB"; cpus params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input: tuple path(fasta_file), val(ext); path hit_ids
    // 修正: 出力ファイル名に拡張子を明示的に含めて `.faa` と `.fna` の競合を防止
    output: tuple val("filtered_catalog"), path("filtered_catalog_${ext}.${ext}")
    script: """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        valid_ids = {line.strip() for line in open("${hit_ids}") if line.strip()}
        with open("${fasta_file}", "r") as fin, open("filtered_catalog_${ext}.${ext}", "w") as fout:
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
    input: path id_mapping; path hit_ids
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
// 1. サブワークフロー
// ==========================================

workflow BACTERIOME_PIPELINE_SUB {
    take: p; host_ref; ref_or_db; taxid_map; ko_map; ko_name_map; taxid_name_map; reads

    main:
    // 1. FASTP による QC・トリミング
    fp = FASTP_SUB(p, reads)
    fp.out.view { "DEBUG [FASTP]: $it" }
    
    // 2. ホスト除去
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)
    host_removed.reads.view { "DEBUG [REMOVE_HOST]: $it" }
    
    // 3. アセンブリ
    asm_res = ASSEMBLY_SUB(p, host_removed.reads)
    asm_res.asm.view { "DEBUG [ASSEMBLY asm]: $it" }
    asm_res.flt_seqs.view { "DEBUG [ASSEMBLY flt_seqs]: $it" }
    
    // 4. Prodigal による ORF (遺伝子) 予測
    orf_res = PRODIGAL_SUB(p, asm_res.flt_seqs)
    orf_res.out.view { "DEBUG [PRODIGAL]: $it" }

    // 5. merge_fasta 用に全サンプルの faa / fna リストを集約
    sample_orf_fna_ch = orf_res.out.map { qry_id, faa, fna, gbk -> tuple(qry_id, fna) }
    sample_orf_fna_ch.view{ "DEBUG [sample_orf_fna_ch]: $it" }
    all_faa_files = orf_res.out.map { qry_id, faa, fna, gbk -> faa }.collect()
    all_faa_files.view{ "DEBUG [all_faa_files]: $it" }
    all_fna_files = orf_res.out.map { qry_id, faa, fna, gbk -> fna }.collect()
    all_fna_files.view{ "DEBUG [all_fna_files]: $it" }

    merge_inputs_ch = all_faa_files.map { tuple(it, "faa") }.mix(all_fna_files.map { tuple(it, "fna") })
    merge_inputs_ch.view{ "DEBUG [merge_inputs_ch]: $it" }
    merged_fastas = merge_fasta(merge_inputs_ch)
    merged_fastas.view{ "DEBUG [merged_fastas]: $it" }

    merged_faa_fasta = merged_fastas.filter { id, path -> path.name.endsWith('.faa') }
    merged_faa_fasta.view{ "DEBUG [merged_faa_fastas]: $it" }
    merged_fna_fasta = merged_fastas.filter { id, path -> path.name.endsWith('.fna') }
    merged_fna_fasta.view{ "DEBUG [merged_fna_fastas]: $it" }

    // 6. NRカタログ構築 (MMseqs2によるクラスタリング ＆ NRCAT連番ID化 ＆ 対応表作成)
    nr_catalog_res = NR_CATALOG_SUB(p, merged_faa_fasta, merged_fna_fasta)
    nr_catalog_res.rep_faa.view { "DEBUG [NR_CATALOG rep_faa]: $it" }
    nr_catalog_res.rep_fna.view { "DEBUG [NR_CATALOG rep_fna]: $it" }
    nr_catalog_res.mapping.view { "DEBUG [NR_CATALOG mapping]: $it" }
    
    nr_rep_faa_path = nr_catalog_res.rep_faa.map { id, path -> path }
    nr_rep_faa_path.view{ "DEBUG [nr_rep_faa_path]: $it" }
    nr_rep_fna_path = nr_catalog_res.rep_fna.map { id, path -> path }
    nr_rep_fna_path.view{ "DEBUG [nr_rep_fna_path]: $it" }
    
    // 7. アノテーション
    annotation_res = ANNOTATE_CATALOG_SUB(p, ref_or_db, nr_catalog_res.rep_faa, taxid_map, ko_map)
    annotation_res.annotated.view { "DEBUG [ANNOTATE_CATALOG annotated]: $it" }
    annotation_res.hit_ids.view { "DEBUG [ANNOTATE_CATALOG hit_ids]: $it" }

    annotated_path_ch = annotation_res.annotated.map { id, path -> path }
    annotated_path_ch.view { "DEBUG [ANNOTATED]: $it" }
    hit_ids_path_ch   = annotation_res.hit_ids.map { id, path -> path }
    hit_ids_path_ch.view { "DEBUG [HIT_IDS]: $it" }
    
    // 8. hit_ids に基づいてカタログ等をフィルタリング
    catalog_to_filter = nr_rep_faa_path.map { tuple(it, "faa") }.mix(nr_rep_fna_path.map { tuple(it, "fna") })
    filtered_catalogs = filter_catalog_by_hits(catalog_to_filter, hit_ids_path_ch)
    filtered_catalogs.view { "DEBUG [FILTERED_CATALOGS]: $it" }

    filtered_rep_faa = filtered_catalogs.filter { id, path -> path.name.endsWith('.faa') }
    filtered_rep_faa.view { "DEBUG [FILTERED_REP_FAA]: $it" }
    filtered_rep_fna = filtered_catalogs.filter { id, path -> path.name.endsWith('.fna') }
    filtered_rep_fna.view { "DEBUG [FILTERED_REP_FNA]: $it" }
    
    mapping_path_ch = nr_catalog_res.mapping.map { id, path -> path }
    filtered_mapping_res = filter_mapping_by_hits(mapping_path_ch, hit_ids_path_ch)
    filtered_mapping_res.view { "DEBUG [FILTERED_MAPPING]: $it" }

    // 9. ORF マッピング
    samples_mapping_res = ALIGN_SAMPLES_TO_CATALOG_SUB(p, filtered_rep_fna, sample_orf_fna_ch)
    samples_mapping_res.out.view { "DEBUG [SAMPLES_MAPPING]: $it" }
    
    // 10. ORF -> TaxID / KO 対応テーブルの構築
    samples_anno_res = ANNOTATE_SAMPLES_SUB(p, samples_mapping_res.out, annotated_path_ch)
    samples_anno_res.samples_anno.view { "DEBUG [SAMPLES_ANNO]: $it" }
    
    // 11. サンプルごとの集計処理 
    summary_res = SUMMARIZE_KO_TAXID_SUB(p, samples_anno_res.samples_anno, ko_name_map, taxid_name_map)
    summary_res.ko_summary.view { "DEBUG [KO_SUMMARY]: $it" }
    summary_res.tax_summary.view { "DEBUG [TAX_SUMMARY]: $it" }

    emit:
    raw_reads              = reads
    fastp_reads            = fp.out
    host_removed_reads     = host_removed.reads
    contigs                = asm_res.asm
    flt_seqs               = asm_res.flt_seqs
    orfs                   = orf_res.out
    nr_catalog_faa         = nr_catalog_res.rep_faa
    nr_catalog_fna         = nr_catalog_res.rep_fna
    annotated              = annotated_path_ch
    hit_ids                = hit_ids_path_ch
    filtered_rep_faa       = filtered_rep_faa
    filtered_rep_fna       = filtered_rep_fna
    filtered_mapping       = filtered_mapping_res
    samples_mapping        = samples_mapping_res.out
    samples_annotation     = samples_anno_res.samples_anno
    ko_summary             = summary_res.ko_summary
    tax_summary            = summary_res.tax_summary    
}

workflow BACTERIOME_PIPELINE_ALL {
    // ==========================================
    // 1. 必須パラメータ・ファイルの厳格バリデーション
    // ==========================================
    if (!params.bacteriome_pipeline_reads) {
        error "【パラメータ不足】 '--bacteriome_pipeline_reads' が指定されていません。"
    }
    if (!params.remove_host_ref_fasta_or_db) {
        error "【パラメータ不足】 '--remove_host_ref_fasta_or_db' が指定されていません。"
    }
    if (!params.annotate_catalog_prot_ref_or_db) {
        error "【パラメータ不足】 '--annotate_catalog_prot_ref_or_db' が指定されていません。"
    }

    def checkMandatoryFile = { paramName, pathStr ->
        if (!pathStr) {
            error "【パラメータ不足】 '--${paramName}' が指定されていません。"
        }
        def f = file(pathStr)
        if (!f.exists()) {
            error "【ファイル非存在エラー】 '--${paramName}' で指定されたファイルが存在しません: ${pathStr}"
        }
        return f
    }

    taxid_map      = checkMandatoryFile('taxid_map_path', params.taxid_map_path)
    ko_map         = checkMandatoryFile('ko_map_path', params.ko_map_path)
    ko_name_map    = checkMandatoryFile('ko_name_map_path', params.ko_name_map_path)
    taxid_name_map = checkMandatoryFile('taxid_name_map_path', params.taxid_name_map_path)

    // ==========================================
    // 2. データフロー構築
    // ==========================================
    p = createNullParamsChannel()

    def reads_list = params.bacteriome_pipeline_reads.split(';')
    def individual_channels = reads_list.collect { Channel.fromFilePairs(it, checkIfExists: true) }
    def reads_mixed = individual_channels.first()
    if (individual_channels.size() > 1) { 
        individual_channels.tail().each { reads_mixed = reads_mixed.mix(it) } 
    }
    
    def index = 0
    reads = reads_mixed.map { id, pair -> tuple("${String.format('%02d', index++)}_${id}", pair) }

    host_ref  = createSeqsChannel(params.remove_host_ref_fasta_or_db)
    ref_or_db = createSeqsChannel(params.annotate_catalog_prot_ref_or_db)

    // サブワークフロー呼び出し
    out_ch = BACTERIOME_PIPELINE_SUB(p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads)
}

workflow { BACTERIOME_PIPELINE_ALL() }