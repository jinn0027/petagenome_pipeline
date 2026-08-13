#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def MERGE_FASTA_MAX_MEMORY  = 16
def MERGE_FASTA_MAX_THREADS = 4

// 3. 上限値による動的クリッピング
params.merge_fasta_memory  = Math.min(params.memory as Integer, MERGE_FASTA_MAX_MEMORY)
params.merge_fasta_threads = Math.min(params.threads as Integer, MERGE_FASTA_MAX_THREADS)

include { createNullParamsChannel; createSeqsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { FASTP_SUB }                 from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB }           from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { ASSEMBLY_SUB }              from "${params.petagenomeDir}/nf/lv2/assembly.nf"
include { PRODIGAL_SUB }              from "${params.petagenomeDir}/nf/lv1/prodigal.nf"
include { ANNOTATE_ORFS_SUB }         from "${params.petagenomeDir}/nf/lv2/annotate_orfs_prot.nf"
include { ALIGN_SAMPLES_TO_ORFS_SUB } from "${params.petagenomeDir}/nf/lv2/align_samples_to_orfs.nf"
include { ANNOTATE_SAMPLES_SUB }      from "${params.petagenomeDir}/nf/lv2/annotate_samples.nf"
include { SUMMARIZE_KO_TAXID_SUB }    from "${params.petagenomeDir}/nf/lv2/summarize_ko_taxid.nf"

// 全サンプルの FASTA (.faa や .fna) を 1 つに結合する汎用プロセス
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
        // stageAs: "?/*" を指定してファイル名の重複衝突を防ぐ
        tuple path(fasta_files, stageAs: "?/*"), val(ext)

    output:
        tuple val("merged_all_samples"), path("combined_orfs.${ext}"), emit: merged_fasta

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    cat ${fasta_files} > combined_orfs.${ext}
    """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

workflow BACTERIOME_PIPELINE_SUB {
    take:
    p
    host_ref       // ホストゲノム(FASTAまたはDB)
    ref_or_db      // アノテーション用リファレンス
    taxid_map      // uniprot_to_taxid.tsv
    ko_map         // uniprot_to_ko.tsv
    ko_name_map    // ko_to_name.tsv (任意)
    taxid_name_map // taxid_to_name.tsv (任意)
    reads          // 入力リードペア [ pair_id, [R1, R2] ]

    main:
    // 1. FASTP による QC・トリミング
    fp = FASTP_SUB(p, reads)

    // 2. ホスト除去
    host_removed = REMOVE_HOST_SUB(p, host_ref, fp.out)

    // 3. アセンブリ
    l_thre = params.containsKey('assembly_l_thre') ? params.assembly_l_thre : 1000
    asm_res = ASSEMBLY_SUB(p, host_removed.reads, l_thre)

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

    // merge_fasta を 1 度呼び出し（内部で並列 2 タスクが起動）
    merged_fastas = merge_fasta(merge_inputs_ch)

    // 出力結果を .faa と .fna に分岐
    merged_faa_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.faa') }
    merged_fna_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.fna') }

    // アノテーション（タンパク質 .faa 側）
    annotation_res = ANNOTATE_ORFS_SUB(p, ref_or_db, merged_faa_fasta, taxid_map, ko_map)

    // ORF マッピング（各サンプルの ORF vs 統合 ORF カタログ）
    samples_mapping_res = ALIGN_SAMPLES_TO_ORFS_SUB(p, merged_fna_fasta, sample_orf_fna_ch)

    // 6. ORF -> TaxID / KO 対応テーブルの構築
    // ANNOTATE_SAMPLES_SUB 内で 1 vs N (samples_mapping_res.out vs annotation_res.annotated) の結合を行う
    samples_anno_res = ANNOTATE_SAMPLES_SUB(
        p,
        samples_mapping_res.out,
        annotation_res.annotated
    )

    // 7. サンプルごとの KO / TaxID 比率の集計処理 
    summary_res = SUMMARIZE_KO_TAXID_SUB(p, samples_anno_res.samples_anno, ko_name_map, taxid_name_map)

    emit:
    raw_reads          = reads
    fastp_reads        = fp.out
    host_removed_reads = host_removed.reads
    contigs            = asm_res.asm
    flt_seqs           = asm_res.flt_seqs
    stats              = asm_res.sts
    orfs               = orf_res.out
    annotated          = annotation_res.annotated
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
    ref_or_db = createSeqsChannel(params.annotate_orfs_prot_ref_or_db)
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
    out_ch.samples_map_out.view    { i -> "SAMPLES MAP RESULT  : $i" }
    out_ch.samples_anno.view       { i -> "SAMPLES KO TAXID TSV: $i" }
    out_ch.ko_summary.view         { i -> "KO SUMMARY TSV      : $i" }
    out_ch.tax_summary.view        { i -> "TAXID SUMMARY TSV   : $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}