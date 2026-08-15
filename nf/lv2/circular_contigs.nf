#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4

// 2. プロセス固有の上限設定
def CLASSIFY_MAX_MEMORY     = 16
def CLASSIFY_MAX_THREADS    = 8

def DEDUPLICATE_MAX_MEMORY  = 8
def DEDUPLICATE_MAX_THREADS = 4

// 3. 上限値による動的クリッピング
params.circular_contigs_classify_memory    = Math.min(params.memory as Integer, CLASSIFY_MAX_MEMORY)
params.circular_contigs_classify_threads   = Math.min(params.threads as Integer, CLASSIFY_MAX_THREADS)

params.circular_contigs_deduplicate_memory = Math.min(params.memory as Integer, DEDUPLICATE_MAX_MEMORY)
params.circular_contigs_deduplicate_threads= Math.min(params.threads as Integer, DEDUPLICATE_MAX_THREADS)

// スイッチ・各種しきい値デフォルト設定
params.use_pzlast           = false
params.pzlast_cfg1          = null
params.pzlast_cfg2          = null

params.circular_contigs_e_thre     = "1e-10"
params.circular_contigs_pi_thre    = "100"
params.circular_contigs_al_thre    = "50"
params.circular_contigs_pi_thre_rd = "95"
params.circular_contigs_qc_thre_rd = "95"
params.circular_contigs_len_l      = "5000"
params.circular_contigs_len_c      = "1500"
params.circular_contigs_pi_self    = 100
params.circular_contigs_al_self    = 50

params.circular_contigs_blast1_num_alignments = 5
params.circular_contigs_blast2_num_alignments = 50

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb as blast_makerefdb1; blastn as blastn1 } from "${params.petagenomeDir}/nf/lv1/blast"
include { blast_makerefdb as blast_makerefdb2; blastn as blastn2 } from "${params.petagenomeDir}/nf/lv1/blast"

if (params.use_pzlast) {
    include { pzlast_makerefdb as pzlast_makerefdb1; pzlastn as pzlastn1 } from "${params.pzrepoDir}/nf/lv1/pzlast"
    include { pzlast_makerefdb as pzlast_makerefdb2; pzlastn as pzlastn2 } from "${params.pzrepoDir}/nf/lv1/pzlast"
}

// ==========================================
// 1. プロセス定義
// ==========================================

process classify {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.circular_contigs_classify_memory}"
    def threads = "${params.circular_contigs_classify_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(ref_id), val(qry_id), path(in_tsv, arity: '1')
        tuple val(qry_id), path(qry, arity: '1')

    output:
        tuple val(qry_id),
              path("${qry_id}/circular.cut.fa"),
              path("${qry_id}/circular.extended.fa"),
              path("${qry_id}/circular.fa"),
              path("${qry_id}/linear.fa"),
              path("${qry_id}/selfhit.tsv", arity: '1')

    script:
    def al_self = getParam(p, params, 'circular_contigs_al_self')
    """
    echo "${processProfile(task)}" | tee prof.txt
    
    qry_=${qry}
    if [[ "\${qry_}" == *.gz ]]; then
        qry_=\${qry_%.gz}
        unpigz -p ${threads} -c ${qry} > \${qry_}
    fi
    mkdir -p ${qry_id}

    # 1. 自己ヒットの抽出
    awk -F "\\t" 'BEGIN{OFS="\\t"} \$1==\$2 {print \$0}' ${in_tsv} > ${qry_id}/selfhit.tsv

    # 2. AWKで一括解析し、各コンティグの最大アライメント末尾座標(pos_end)を抽出
    awk -F "\\t" -v al_self="${al_self}" '
    BEGIN { OFS="\\t" }
    {
        q1 = (\$7 < \$8) ? \$7 : \$8;
        q2 = (\$7 < \$8) ? \$8 : \$7;
        t1 = (\$9 < \$10) ? \$9 : \$10;
        t2 = (\$9 < \$10) ? \$10 : \$9;

        # 非自己アライメント領域のチェック
        if (q1 != t1 || q2 != t2) {
            aln_len = \$4;
            if (aln_len >= al_self) {
                pos = 0;
                if (q1 == 1) pos = t1;
                else if (t1 == 1) pos = q1;

                if (pos > 1) {
                    if (pos > max_pos[\$1]) max_pos[\$1] = pos;
                }
            }
        }
    }
    END {
        for (id in max_pos) {
            print id, max_pos[id];
        }
    }' ${qry_id}/selfhit.tsv > ${qry_id}/circular_positions.txt

    # 3. ID リストと BED ファイルの作成
    cut -f1 ${qry_id}/circular_positions.txt > ${qry_id}/circular_ids.txt
    awk 'BEGIN{OFS="\\t"} {print \$1, 1, \$2}' ${qry_id}/circular_positions.txt > ${qry_id}/cut.bed

    # 4. SeqKit による一括分画と一括切断処理 (Bashループなし)
    if [ -s ${qry_id}/circular_ids.txt ]; then
        seqkit grep -f ${qry_id}/circular_ids.txt \${qry_} > ${qry_id}/circular.fa
        seqkit grep -v -f ${qry_id}/circular_ids.txt \${qry_} > ${qry_id}/linear.fa

        # 切断処理（seqkit seq -i でヘッダー末尾の座標テキストをカットしIDを保持）
        seqkit subseq --bed ${qry_id}/cut.bed ${qry_id}/circular.fa | seqkit seq -i > ${qry_id}/circular.cut.fa
        cat ${qry_id}/circular.cut.fa ${qry_id}/circular.fa > ${qry_id}/circular.extended.fa
    else
        touch ${qry_id}/circular.fa ${qry_id}/circular.cut.fa ${qry_id}/circular.extended.fa
        cp \${qry_} ${qry_id}/linear.fa
    fi
    """
}

process deduplicate {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }

    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb      = "${params.circular_contigs_deduplicate_memory}"
    def threads = "${params.circular_contigs_deduplicate_threads}"

    memory params.executor == "sge" ? null : "${gb} GB"
    cpus   params.executor == "sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
        tuple val(p), val(id), path(circular_cut), path(circular_ext), path(circular), path(blst_out_tsv, arity: '1')

    output:
        tuple val(id),
              path("${id}/circular.cut.fa"),
              path("${id}/circular.extended.fa"),
              path("${id}/circular.fa"),
              path("${id}/otherhit.tsv", arity: '1'),
              path("${id}/*.txt", arity: '3')

    script:
    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${id}
    awk -F "\\t" '{OFS="\\t"} { if (\$1 != \$2) print \$0 }' ${blst_out_tsv} > ${id}/otherhit.tsv
    python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py ${circular_cut} > ${id}/circular_cut.all.length.txt
    ruby ${params.petagenomeDir}/scripts/Ruby/extract_contig_redundancy.3.rb -b ${id}/otherhit.tsv -l ${id}/circular_cut.all.length.txt -c ${params.circular_contigs_qc_thre_rd} -d 6 -i ${id}/out_rd_info.txt -o ${id}/out_ex_config.txt
    touch ${id}/circular.cut.fa ${id}/circular.extended.fa ${id}/circular.fa
    python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_cut} | sed '/^\$/d' > ${id}/circular.cut.fa
    python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_ext} | sed '/^\$/d' > ${id}/circular.extended.fa
    python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular} | sed '/^\$/d' > ${id}/circular.fa
    """
}

// ==========================================
// 2. サブワークフロー（本体）
// ==========================================

workflow CIRCULAR_CONTIGS_SUB {
    take:
    p
    contig

    main:
    
    def use_pzlast_closure = { getParam(p, params, 'use_pzlast') }
    if (use_pzlast_closure()) {
        // --- PZLAST モード ---
        cfg1 = params.pzlast_cfg1 ? Channel.of(file(params.pzlast_cfg1)) : Channel.of('')

        // 1. PZLAST DB1 作成 ＋ アライメント 1
        db1_in = p.combine(contig).combine(cfg1).map { p_val, id, seq, cfg ->
            tuple(p_val, id, seq, cfg)
        }
        pzlstdb1 = pzlast_makerefdb1(db1_in)

        pzlstin1 = pzlstdb1.combine(contig).map { ref_id, db_path, qry_id, seq_path ->
            tuple(ref_id, db_path, qry_id, seq_path)
        }
        p_pzlastn1 = Channel.of([
            'pzlast_outfmt': 6,
            'pzlast_fmt6_swapside': 's',
            'pzlast_q_with_comp': 0
        ])
        
        aln1_in = p_pzlastn1.combine(pzlstin1).map { p_map, ref_id, db_path, qry_id, seq_path ->
            tuple(p_map, ref_id, db_path, qry_id, seq_path)
        }
        aln1 = pzlastn1(aln1_in).combine(cfg1)

        // 2. 分類 (Classify)
        clsfy_p_in = p.combine(aln1).map { p_val, ref_id, qry_id, tsv, cfg ->
            tuple(p_val, ref_id, qry_id, tsv)
        }
        clsfy = classify(clsfy_p_in, contig)

        circular_cut = clsfy.map { qry_id, cut, ext, circular, linear, selfhit_tsv ->
            tuple(qry_id, cut)
        }

        // 3. PZLAST DB2 作成 ＋ アライメント 2
        cfg2 = params.pzlast_cfg2 ? Channel.of(file(params.pzlast_cfg2)) : Channel.of('')

        db2_in = p.combine(circular_cut).combine(cfg2).map { p_val, id, seq, cfg ->
            tuple(p_val, id, seq, cfg)
        }
        pzlstdb2 = pzlast_makerefdb2(db2_in)

        p_pzlastn2 = Channel.of([
            'pzlast_outfmt': 6,
            'pzlast_fmt6_swapside': 's',
            'pzlast_q_with_comp': 1
        ])
        pzlstin2 = pzlstdb2.combine(circular_cut).map { ref_id, db_path, qry_id, seq_path ->
            tuple(ref_id, db_path, qry_id, seq_path)
        }
        aln2_in = p_pzlastn2.combine(pzlstin2).map { p_map, ref_id, db_path, qry_id, seq_path ->
            tuple(p_map, ref_id, db_path, qry_id, seq_path)
        }
        aln2 = pzlastn2(aln2_in).combine(cfg2)

        // 4. 重複除去 (Deduplicate)
        // aln2 の要素数が異なる場合にも対応できるようインデックスで参照
        aln2_key  = aln2.map { item -> tuple(item[1], item[2]) } // [qry_id, tsv]
        clsfy_key = clsfy.map { qry_id, cut, ext, circular, linear, selfhit_tsv -> tuple(qry_id, cut, ext, circular) }

        ch_new = clsfy_key.join(aln2_key).map { qry_id, cut, ext, circular, aln2_tsv ->
            tuple(qry_id, cut, ext, circular, aln2_tsv)
        }

        dedupl_in = p.combine(ch_new).map { p_val, qry_id, cut, ext, circular, aln2_tsv ->
            tuple(p_val, qry_id, cut, ext, circular, aln2_tsv)
        }
        dedupl = deduplicate(dedupl_in)

    } else {
        // --- BLAST モード ---
        // 1. BLAST DB1 作成 ＋ アライメント 1
        db1_in = p.combine(contig).map { p_val, id, seq ->
            tuple(p_val, id, seq)
        }
        blstdb1 = blast_makerefdb1(db1_in)

        blstin1 = blstdb1.combine(contig).map { ref_id, db_path, qry_id, seq_path ->
            tuple(ref_id, db_path, qry_id, seq_path)
        }
        p_blastn1 = Channel.of([
            'blast_task': 'megablast',
            'blast_perc_identity': params.circular_contigs_pi_self,
            'blast_evalue': params.circular_contigs_e_thre,
            'blast_outfmt': 6,
            'blast_num_alignments': params.circular_contigs_blast1_num_alignments,
            'blast_strand': 'plus'
        ])
        
        aln1_in = p_blastn1.combine(blstin1).map { p_map, ref_id, db_path, qry_id, seq_path ->
            tuple(p_map, ref_id, db_path, qry_id, seq_path)
        }
        aln1 = blastn1(aln1_in)

        // 2. 分類 (Classify)
        clsfy_p_in = p.combine(aln1).map { p_val, ref_id, qry_id, tsv ->
            tuple(p_val, ref_id, qry_id, tsv)
        }
        clsfy = classify(clsfy_p_in, contig)

        circular_cut = clsfy.map { qry_id, cut, ext, circular, linear, selfhit_tsv ->
            tuple(qry_id, cut)
        }

        // 3. BLAST DB2 作成 ＋ アライメント 2
        db2_in = p.combine(circular_cut).map { p_val, id, seq ->
            tuple(p_val, id, seq)
        }
        blstdb2 = blast_makerefdb2(db2_in)

        p_blastn2 = Channel.of([
            'blast_task': 'megablast',
            'blast_perc_identity': params.circular_contigs_pi_self,
            'blast_evalue': params.circular_contigs_e_thre,
            'blast_outfmt': 6,
            'blast_num_alignments': params.circular_contigs_blast2_num_alignments,
            'blast_strand': 'both'
        ])
        blstin2 = blstdb2.combine(circular_cut).map { ref_id, db_path, qry_id, seq_path ->
            tuple(ref_id, db_path, qry_id, seq_path)
        }
        aln2_in = p_blastn2.combine(blstin2).map { p_map, ref_id, db_path, qry_id, seq_path ->
            tuple(p_map, ref_id, db_path, qry_id, seq_path)
        }
        aln2 = blastn2(aln2_in)

        // 4. 重複除去 (Deduplicate)
        aln2_key  = aln2.map { item -> tuple(item[1], item[2]) } // [qry_id, tsv]
        clsfy_key = clsfy.map { qry_id, cut, ext, circular, linear, selfhit_tsv -> tuple(qry_id, cut, ext, circular) }

        ch_new = clsfy_key.join(aln2_key).map { qry_id, cut, ext, circular, blst2_tsv ->
            tuple(qry_id, cut, ext, circular, blst2_tsv)
        }

        dedupl_in = p.combine(ch_new).map { p_val, qry_id, cut, ext, circular, aln2_tsv ->
            tuple(p_val, qry_id, cut, ext, circular, aln2_tsv)
        }
        dedupl = deduplicate(dedupl_in)
    }

    emit:
    aln1   = aln1
    clsfy  = clsfy
    aln2   = aln2
    dedupl = dedupl
}

// ==========================================
// 3. コマンドライン (-entry) 用エントリーポイント
// ==========================================

workflow CIRCULAR_CONTIGS_ALL {
    p      = createNullParamsChannel()
    contig = createSeqsChannel(params.circular_contigs_contig)

    out_ch = CIRCULAR_CONTIGS_SUB(p, contig)

    out_ch.aln1.view   { i -> "ALN1: $i" }
    out_ch.clsfy.view  { i -> "CLSFY: $i" }
    out_ch.aln2.view   { i -> "ALN2: $i" }
    out_ch.dedupl.view { i -> "DEDUPL: $i" }
}

workflow {
    CIRCULAR_CONTIGS_ALL()
}