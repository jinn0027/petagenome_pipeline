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
include { SUMMARIZE_KO_TAXID_SUB }    from "${params.petagenomeDir}/nf/lv2/summarize_ko_taxid.nf"

// 全サンプルの FASTA (.faa や .fna) を 1 つに結合する汎用プロセス
process MERGE_FASTA {
    tag "${ext}"

    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions}"
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    input:
    tuple path(fasta_files), val(ext) // (配列ファイル群, 拡張子) のペアとして受け取る

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

    // --------------------------------------------------------------------------
    // 各サンプル固有の ORF 塩基配列 (.fna) チャンネルを作成
    // orf_res.out: tuple(qry_id, faa_path, fna_path, gbk_path)
    // --------------------------------------------------------------------------
    sample_orf_fna_ch = orf_res.out.map { qry_id, faa, fna, gbk -> tuple(qry_id, fna) }

    // 5. PRODIGALのサンプル毎の結果をまとめる

    // .faa ファイルのリストを作成（この時、qry_id を使ってユニークな名前へリネーム）
    all_faa_files = orf_res.out
        .map { qry_id, faa, fna, gbk -> 
            return faa.moveTo("${faa.parent}/${qry_id}.faa")
        }
        .collect()

    // .fna ファイルのリストを作成（同様にユニーク名へリネーム）
    all_fna_files = orf_res.out
        .map { qry_id, faa, fna, gbk -> 
            return fna.moveTo("${fna.parent}/${qry_id}.fna")
        }
        .collect()

    // MERGE_FASTA への入力チャンネル作成
    merge_inputs_ch = all_faa_files.map { faa_list -> tuple(faa_list, "faa") }
        .mix(all_fna_files.map { fna_list -> tuple(fna_list, "fna") })
    // MERGE_FASTA を 1 度呼び出し（内部で並列 2 タスクが起動）
    merged_fastas = MERGE_FASTA(merge_inputs_ch)

    // 出力結果を .faa と .fna に分岐
    merged_faa_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.faa') }
    merged_fna_fasta = merged_fastas.merged_fasta.filter { id, path -> path.name.endsWith('.fna') }

    // 6. Contig/ORF -> TaxID / KO 対応テーブルの構築

    // アノテーション（タンパク質 .faa 側）
    annotation_res = ANNOTATE_TAXID_KO_SUB(p, ref_or_db, merged_faa_fasta, taxid_map, ko_map)

    // クエリにsample_orf_fna_ch (ORF) を渡す
    contig_mapping_res = ALIGN_CONTIGS_TO_ORFS_SUB(p, merged_fna_fasta, sample_orf_fna_ch)

    contig_anno_res = ANNOTATE_CONTIG_SUB(p, contig_mapping_res.out, annotation_res.annotated)

    // 7. サンプルごとの KO / TaxID 比率の集計処理
    summary_res = SUMMARIZE_KO_TAXID_SUB(p, contig_anno_res.contig_anno)

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
    out_ch.ko_summary.view         { i -> "KO SUMMARY TSV     : $i" }
    out_ch.tax_summary.view        { i -> "TAXID SUMMARY TSV  : $i" }
}

// デフォルトエントリーポイント
workflow {
    BACTERIOME_PIPELINE_ALL()
}