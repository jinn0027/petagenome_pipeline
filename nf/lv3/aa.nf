#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def MERGE_FASTA_MAX_MEMORY  = 16
def MERGE_FASTA_MAX_THREADS = 4
def FILTER_CATALOG_MAX_MEMORY  = 8
def FILTER_CATALOG_MAX_THREADS = 2

// 3. 上限値による動的クリッピング
params.merge_fasta_memory   = Math.min(params.memory as Integer, MERGE_FASTA_MAX_MEMORY)
params.merge_fasta_threads  = Math.min(params.threads as Integer, MERGE_FASTA_MAX_THREADS)
params.filter_catalog_memory   = Math.min(params.memory as Integer, FILTER_CATALOG_MAX_MEMORY)
params.filter_catalog_threads  = Math.min(params.threads as Integer, FILTER_CATALOG_MAX_THREADS)

include { createNullParamsChannel; createSeqsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { FASTP_SUB }                  from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB }              from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { ASSEMBLY_SUB }                 from "${params.petagenomeDir}/nf/lv2/assembly.nf"
include { PRODIGAL_SUB }                 from "${params.petagenomeDir}/nf/lv1/prodigal.nf"
include { ANNOTATE_CATALOG_SUB }         from "${params.petagenomeDir}/nf/lv2/annotate_catalog_prot.nf"
include { ALIGN_SAMPLES_TO_CATALOG_SUB } from "${params.petagenomeDir}/nf/lv2/align_samples_to_catalog.nf"
include { ANNOTATE_SAMPLES_SUB }         from "${params.petagenomeDir}/nf/lv2/annotate_samples.nf"
include { SUMMARIZE_KO_TAXID_SUB }       from "${params.petagenomeDir}/nf/lv2/summarize_ko_taxid.nf"

// 新しく作成したカタログ構築サブワークフローをインポート
include { NR_CATALOG_SUB }          from "${params.petagenomeDir}/nf/lv2/nr_catalog.nf"

// 全サンプルの FASTA (.faa や .fna) を単純に 1 つに結合する汎用プロセス（カタログ作成前の入力用）
process merge_fasta {
    tag "${ext}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.merge_fasta_memory}"
    def threads = "${params.merge_fasta_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple path(fasta_files, stageAs: "?/*"), val(ext)

    output:
        tuple val("merged_all_samples"), path("combined_orfs.${ext}"), emit: merged_fasta

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        cat ${fasta_files} > combined_orfs.${ext}
        """
}

// hit_ids に含まれる ID のみにカタログ（FAA/FNA）をフィルタリングするプロセス
process filter_catalog_by_hits {
    tag "${ext}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.filter_catalog_memory}"
    def threads = "${params.filter_catalog_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple path(fasta_file), val(ext)
        path hit_ids

    output:
        tuple val("filtered_catalog"), path("filtered_catalog.${ext}"), emit: filtered_fasta

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        # 1. hit_ids の読み込み
        valid_ids = set()
        with open("${hit_ids}", "r") as f:
            for line in f:
                qid = line.strip()
                if qid:
                    valid_ids.add(qid)

        # 2. FASTAのフィルタリング
        input_fasta = "${fasta_file}"
        output_fasta = "filtered_catalog.${ext}"

        with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
            write_flag = False
            for line in fin:
                if line.startswith(">"):
                    header = line.strip()[1:].split()[0]
                    if header in valid_ids:
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

// アノテーションの hit_ids に基づいて id_mapping.tsv をフィルタリングするプロセス
process filter_mapping_by_hits {
    tag "filter_mapping"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.filter_catalog_memory}"
    def threads = "${params.filter_catalog_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        path id_mapping
        path hit_ids

    output:
        tuple val("filtered_mapping"), path("filtered_id_mapping.tsv"), emit: filtered_mapping

    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        python3 - << 'EOF'
        # 1. hit_ids の読み込み
        valid_ids = set()
        with open("${hit_ids}", "r") as f:
            for line in f:
                qid = line.strip()
                if qid:
                    valid_ids.add(qid)

        # 2. id_mapping.tsv のフィルタリング (形式: new_id \t old_id)
        with open("${id_mapping}", "r") as fin, open("filtered_id_mapping.tsv", "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if not parts:
                    continue
                catalog_id = parts[0]
                if catalog_id in valid_ids:
                    fout.write(line)
        EOF
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

workflow BACTERIOME_PIPELINE_SUB {
    take:
    p
    host_ref           // ホストゲノム(FASTAまたはDB)
    ref_or_db          // アノテーション用リファレンス
    taxid_map          // uniprot_to_taxid.tsv
    ko_map             // uniprot_to_ko.tsv
    ko_name_map        // ko_to_name.tsv (任意)
    taxid_name_map     // taxid_to_name.tsv (任意)
    reads              // 入力リードペア [ pair_id, [R1, R2] ]

    main:
    // 1. FASTP による QC・トリミング
    fp = FASTP_SUB(p, reads)

    // 2. ホスト除去
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)

    // 3. アセンブリ
    asm_res = ASSEMBLY_SUB(p, host_removed.reads)

    // 4. Prodigal による ORF (遺伝子) 予測
    orf_res = PRODIGAL_SUB(p, asm_res.flt_seqs)

    // 各サンプル固有の ORF 塩基配列 (.fna) チャンネル（クエリ用）
    sample_orf_fna_ch = orf_res.out.map { qry_id, faa, fna, gbk -> tuple(qry_id, fna) }

    // 5. merge_fasta 用に全サンプルの faa / fna リストを集約
    all_faa_files = orf_res.out.map { qry_id, faa, fna, gbk -> faa }.collect()
    all_fna_files = orf_res.out.map { qry_id, faa, fna, gbk -> fna }.collect()

    // merge_fasta への入力チャンネル作成
    merge_inputs_ch = all_faa_files.map { faa_list -> tuple(faa_list, "faa") }
        .mix(all_fna_files.map { fna_list -> tuple(fna_list, "fna") })

    // merge_fasta を呼び出し（全サンプルのfaa / fnaをそれぞれ結合）
    merged_fastas = merge_fasta(merge_inputs_ch)

    // 出力結果を .faa と .fna に分岐
    merged_faa_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.faa') }
    merged_fna_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.fna') }

    // 6. NRカタログ構築 (MMseqs2によるクラスタリング ＆ 一意なローカルID化 ＆ 対応表作成)
    nr_catalog_res = NR_CATALOG_SUB(p, merged_faa_fasta, merged_fna_fasta)

    // 7. アノテーション（クラスタリング済みの代表アミノ酸カタログ .faa 側に対して実行）
    annotation_res = ANNOTATE_CATALOG_SUB(p, ref_or_db, nr_catalog_res.rep_faa, taxid_map, ko_map)

    // 8. hit_ids に基づいてカタログ（FAAおよびFNA）、マッピングテーブル (id_mapping)をフィルタリング
    catalog_to_filter = nr_catalog_res.rep_faa.map { id, path -> tuple(path, "faa") }
        .mix(nr_catalog_res.rep_fna.map { id, path -> tuple(path, "fna") })

    filtered_catalogs = filter_catalog_by_hits(catalog_to_filter, annotation_res.hit_ids)

    filtered_rep_faa = filtered_catalogs.filtered_fasta.filter { id, path -> path.name.endsWith('.faa') }
    filtered_rep_fna = filtered_catalogs.filtered_fasta.filter { id, path -> path.name.endsWith('.fna') }
    filtered_mapping_res = filter_mapping_by_hits(nr_catalog_res.mapping.map { id, path -> path }, annotation_res.hit_ids)

    // 9. ORF マッピング（各サンプルの ORF vs フィルタリング済みカタログ塩基配列 .fna）
    samples_mapping_res = ALIGN_SAMPLES_TO_CATALOG_SUB(p, filtered_rep_fna, sample_orf_fna_ch)

    // 10. ORF -> TaxID / KO 対応テーブルの構築
    samples_anno_res = ANNOTATE_SAMPLES_SUB(p, samples_mapping_res.out, annotation_res.annotated)

    // 11. サンプルごとの KO / TaxID 比率の集計処理 
    summary_res = SUMMARIZE_KO_TAXID_SUB(p, samples_anno_res.samples_anno, ko_name_map, taxid_name_map)

    emit:
    raw_reads          = reads
    fastp_reads        = fp.out
    host_removed_reads = host_removed.reads
    contigs            = asm_res.asm
    flt_seqs           = asm_res.flt_seqs
    orfs               = orf_res.out
    annotated          = annotation_res.annotated
    hit_ids            = annotation_res.hit_ids
    nr_catalog_faa     = filtered_rep_faa // フィルタリング済みの代表アミノ酸カタログ
    nr_catalog_fna     = filtered_rep_fna // フィルタリング済みの代表塩基カタログ
    id_mapping         = filtered_mapping_res.filtered_mapping // フィルタリング済みのマッピングテーブル
    samples_map_out    = samples_mapping_res.out
    samples_anno       = samples_anno_res.samples_anno
    ko_summary         = summary_res.ko_summary
    tax_summary        = summary_res.tax_summary
}

// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow BACTERIOME_PIPELINE_ALL {
    p = createNullParamsChannel()

    // 1. 入力リードのチャンネル生成
    def reads_list = params.bacteriome_pipeline_reads.split(';')

    def individual_channels = reads_list.collect { reads_path ->
        Channel.fromFilePairs(reads_path, checkIfExists: true)
    }

    def reads_mixed = individual_channels.first()
    if (individual_channels.size() > 1) {
        individual_channels.tail().each { ch ->
            reads_mixed = reads_mixed.mix(ch)
        }
    }

    def index = 0
    reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index++)}_${id}"
        tuple(new_id, pair)
    }

    // 2. 各種リファレンス・マップファイルの準備
    host_ref  = createSeqsChannel(params.remove_host_ref_fasta_or_db)
    ref_or_db = createSeqsChannel(params.annotate_catalog_prot_ref_or_db)
    taxid_map = file(params.taxid_map_path, checkIfExists: true)
    ko_map    = file(params.ko_map_path, checkIfExists: true)

    // 名称マッピングファイル (未指定の場合は NO_FILE ダミーを渡す)
    ko_name_map    = (params.containsKey('ko_name_map_path') && params.ko_name_map_path) ? file(params.ko_name_map_path, checkIfExists: true) : file('NO_FILE')
    taxid_name_map = (params.containsKey('taxid_name_map_path') && params.taxid_name_map_path) ? file(params.taxid_name_map_path, checkIfExists: true) : file('NO_FILE')

    // 3. パイプライン本体の呼び出し
    out_ch = BACTERIOME_PIPELINE_SUB(p, host_ref, ref_or_db, taxid_map, ko_map, ko_name_map, taxid_name_map, reads)

    // 結果の確認
    out_ch.fastp_reads.view        { i -> "FASTP PASSED READS  : $i" }
    out_ch.host_removed_reads.view { i -> "HOST REMOVED READS  : $i" }
    out_ch.contigs.view            { i -> "ASSEMBLY CONTIGS    : $i" }
    out_ch.orfs.view               { i -> "PRODIGAL ORFS       : $i" }
    out_ch.annotated.view          { i -> "ANNOTATED TSV       : $i" }
    out_ch.nr_catalog_faa.view     { i -> "NR CATALOG FAA      : $i" }
    out_ch.nr_catalog_fna.view     { i -> "NR CATALOG FNA      : $i" }
    out_ch.id_mapping.view         { i -> "ID MAPPING TSV      : $i" }
    out_ch.samples_map_out.view    { i -> "SAMPLES MAP RESULT  : $i" }
    out_ch.samples_anno.view       { i -> "SAMPLES KO TAXID TSV: $i" }
    out_ch.ko_summary.view         { i -> "KO SUMMARY TSV      : $i" }
    out_ch.tax_summary.view        { i -> "TAXID SUMMARY TSV   : $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}