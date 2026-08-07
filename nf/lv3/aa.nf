#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { createNullParamsChannel; createSeqsChannel; getParam } from "${params.petagenomeDir}/nf/common/utils"
include { FASTP_SUB }                 from "${params.petagenomeDir}/nf/lv1/fastp.nf"
include { REMOVE_HOST_SUB }           from "${params.petagenomeDir}/nf/lv2/remove_host.nf"
include { ASSEMBLY_SUB }              from "${params.petagenomeDir}/nf/lv2/assembly.nf"
include { PRODIGAL_SUB }              from "${params.petagenomeDir}/nf/lv1/prodigal.nf"
include { ANNOTATE_TAXID_KO_SUB }     from "${params.petagenomeDir}/nf/lv2/annotate_p.nf"
include { ALIGN_CONTIGS_TO_ORFS_SUB } from "${params.petagenomeDir}/nf/lv2/align_contigs_to_orfs.nf"
include { ANNOTATE_CONTIG_SUB }       from "${params.petagenomeDir}/nf/lv2/annotate_contig.nf"

// 全サンプルの FASTA (.faa や .fna) を 1 つに結合する汎用プロセス
process MERGE_FASTA {
    tag "${ext}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output

    input:
    path fasta_files // collect() から渡される全サンプルの FASTA リスト
    val ext          // "faa" または "fna"

    output:
    tuple val("merged_all_samples"), path("combined_orfs.${ext}"), emit: merged_fasta

    script:
    """
    cat ${fasta_files} > combined_orfs.${ext}
    """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

workflow BACTERIOME_PIPELINE_SUB {
    take:
    p
    host_ref      // ホストゲノム(FASTAまたはDB)
    ref_or_db     // アノテーション用リファレンス
    taxid_map     // uniprot_to_taxid.tsv
    ko_map        // uniprot_to_ko.tsv
    reads         // 入力リードペア [ pair_id, [R1, R2] ]

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

    // 5. 全サンプルの .faa (タンパク質) ファイルを集約して結合・アノテーション
    all_faa_files = orf_res.out.map { qry_id, faa, fna, gbk -> faa }.collect()
    merged_faa_res = MERGE_FASTA(all_faa_files, "faa")

    annotation_res = ANNOTATE_TAXID_KO_SUB(p, ref_or_db, merged_faa_res.merged_fasta, taxid_map, ko_map)

    // 6. 全サンプルの .fna (塩基配列) ファイルを集約して結合・Contig マッピング
    all_fna_files = orf_res.out.map { qry_id, faa, fna, gbk -> fna }.collect()
    merged_fna_res = MERGE_FASTA(all_fna_files, "fna")

    // マージした ORF FASTA を DB にビルドし、各サンプルの contig を検索 (PZLAST / MMseqs2)
    contig_mapping_res = ALIGN_CONTIGS_TO_ORFS_SUB(p, merged_fna_res.merged_fasta, asm_res.flt_seqs)

    // 7. Contig -> TaxID / KO 対応テーブルの構築
    contig_anno_res = ANNOTATE_CONTIG_SUB(p, contig_mapping_res.out, annotation_res.annotated)

    emit:
    raw_reads          = reads
    fastp_reads        = fp.out
    host_removed_reads = host_removed.reads
    contigs            = asm_res.asm
    flt_seqs           = asm_res.flt_seqs
    stats              = asm_res.sts
    orfs               = orf_res.out
    annotated          = annotation_res.annotated
    contig_map_out     = contig_mapping_res.out
    contig_taxid_ko    = contig_anno_res.contig_anno
}


// ==========================================
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow BACTERIOME_PIPELINE_ALL {
    p = createNullParamsChannel()

    // 1. 入力リードのチャンネル生成
    def reads_list = params.bacteriome_pipeline_reads.split(';')

    def individual_channels = reads_list.collect { reads_path ->
        channel.fromFilePairs(reads_path, checkIfExists: true)
    }

    def reads_mixed = individual_channels.first()
    individual_channels.tail().each { ch ->
        reads_mixed = reads_mixed.mix(ch)
    }

    def index = 0
    reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index++)}_${id}"
        tuple(new_id, pair)
    }

    // 2. 各種リファレンス・マップファイルの準備
    host_ref  = createSeqsChannel(params.remove_host_ref_fasta_or_db)
    ref_or_db = createSeqsChannel(params.annotate_p_ref_or_db)
    taxid_map = file(params.taxid_map_path, checkIfExists: true)
    ko_map    = file(params.ko_map_path, checkIfExists: true)

    // 3. パイプライン本体の呼び出し
    out_ch = BACTERIOME_PIPELINE_SUB(p, host_ref, ref_or_db, taxid_map, ko_map, reads)

    // 結果の確認
    out_ch.fastp_reads.view        { i -> "FASTP PASSED READS : $i" }
    out_ch.host_removed_reads.view { i -> "HOST REMOVED READS : $i" }
    out_ch.contigs.view            { i -> "ASSEMBLY CONTIGS   : $i" }
    out_ch.orfs.view               { i -> "PRODIGAL ORFS      : $i" }
    out_ch.annotated.view          { i -> "ANNOTATED TSV      : $i" }
    out_ch.contig_map_out.view     { i -> "CONTIG MAP RESULT  : $i" }
    out_ch.contig_taxid_ko.view    { i -> "CONTIG TAXID KO TSV: $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}