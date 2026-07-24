#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.circular_contigs_classify_memory = params.memory
params.circular_contigs_classify_threads = params.threads
params.circular_contigs_deduplicate_memory = params.memory
params.circular_contigs_deduplicate_threads = params.threads

params.circular_contigs_e_thre = "1e-10"        // E-value cutoff for circular formation
params.circular_contigs_pi_thre = "100"         // identity for circular formation
params.circular_contigs_al_thre = "50"          // minimum alignment length for circular formation
params.circular_contigs_pi_thre_rd = "95"       // identity for removal of redundancy
params.circular_contigs_qc_thre_rd = "95"       // query coverage for removal of redundancy
params.circular_contigs_len_l = "5000"          // minimum length for linear contigs
params.circular_contigs_len_c = "1500"          // minimum length for circular contigs
params.circular_contigs_pi_self = 100           // identity for circular formation
params.circular_contigs_al_self = 50            // alignment length for circular formation
params.circular_contigs_blast1_num_alignments=5
params.circular_contigs_blast2_num_alignments=50

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createSeqsChannel } \
    from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb as blast_makerefdb1; blastn as blastn1} from "${params.petagenomeDir}/nf/lv1/blast"
include { blast_makerefdb as blast_makerefdb2; blastn as blastn2} from "${params.petagenomeDir}/nf/lv1/blast"

if (params.use_pzlast) {
  include { pzlast_makerefdb as pzlast_makerefdb1; pzlastn as pzlastn1} from "${params.pzrepoDir}/nf/lv1/pzlast"
  include { pzlast_makerefdb as pzlast_makerefdb2; pzlastn as pzlastn2} from "${params.pzrepoDir}/nf/lv1/pzlast"
}

process classify {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.circular_contigs_classify_memory}"
    def threads = "${params.circular_contigs_classify_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        tuple val(p), val(ref_id), val(qry_id), path(in_tsv, arity: '1')
        tuple val(qry_id), path(qry, arity: '1')
    output:
        tuple val(ref_id),
        path("${qry_id}/circular.cut.fa"),
        path("${qry_id}/circular.extended.fa"),
        path("${qry_id}/circular.fa"),
        path("${qry_id}/linear.fa"),
        path("${qry_id}/selfhit.tsv", arity: '1')
    script:
        """
        echo "${processProfile(task)}" | tee prof.txt
        qry_=${qry}
        echo ${qry} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            qry_=\${qry_%%.gz}
            unpigz -c ${qry} > \${qry_}
        fi
        mkdir -p ${qry_id}

        # 自己ヒットしたもののみを選ぶ
        awk -F "\t" '{OFS="\t"}  {
             if (\$1==\$2) {
                 print \$0
             }
        }' ${in_tsv} > ${qry_id}/selfhit.tsv 

        # 完全に一致する断片（自己アライメント）を除外し、
        # 残ったヒット（環状オーバーラップまたはリピート）を抽出
        awk -F "\t" '{OFS="\t"} {
            q1 = \$7 < \$8 ? \$7 : \$8;
            q2 = \$7 < \$8 ? \$8 : \$7;
            t1 = \$9 < \$(10) ? \$9 : \$(10);
            t2 = \$9 < \$(10) ? \$(10) : \$9;
            if (q1 != t1 || q2 != t2) {
                print \$0
            }
        }' ${qry_id}/selfhit.tsv > ${qry_id}/non_self_aligned_hits.tsv

        # 配列IDごとにヒットを分割
        awk -F "\t" -v qry_id="${qry_id}" '{OFS="\t"} {
            f=\$1; gsub("/", "__", f);
            print \$0 >> qry_id"/selfhit"f".tsv"
        }' ${qry_id}/non_self_aligned_hits.tsv

        touch ${qry_id}/circular.cut.fa ${qry_id}/circular.extended.fa ${qry_id}/circular.fa ${qry_id}/linear.fa

        # 配列IDごとにFASTAファイル分割
        seqkit split -i \${qry_} --by-id-prefix @ -O ${qry_id}

        seqkit fx2tab -j ${params.circular_contigs_classify_threads} -n -i -l \${qry_} | while read -r id len; do
            f="\$(echo \${id} | sed 's#/#__#g')"
            each_fa="${qry_id}/@\${f}.fa"
            each_selfhit="${qry_id}/selfhit"\${f}".tsv"
            if [ -f \${each_selfhit} ] ; then
                pos_end=\$(sort -n -r -k4 \${each_selfhit} | sed -n '1p' | \\
                           awk -v id=\${id} -v len=\${len} -v al_self=${getParam(p, 'circular_contigs_al_self')} '{ \\
                               if (\$1==id && \$4!=len && \$4>=al_self) { \\
                                   q1 = \$7 < \$8 ? \$7 : \$8; \\
                                   t1 = \$9 < \$(10) ? \$9 : \$(10); \\
                                   if (q1==1) {print(t1)} \\
                                   else if (t1==1) {print(q1)}; \\
                               } \\
                           }')
                if [ "\${pos_end}" != "" ] && [ "\${pos_end}" -gt 1 ] ; then
                    cat \${each_fa} | seqkit subseq -r 1:\${pos_end} > _fa
                    cat _fa >> ${qry_id}/circular.cut.fa
                    tail -n +2 \${each_fa} >> _fa
                    cat _fa >> ${qry_id}/circular.extended.fa
                    rm -f _fa
                    cat \${each_fa} >> ${qry_id}/circular.fa
                else
                    cat \${each_fa} >> ${qry_id}/linear.fa
                fi
            fi
        done
        """
}

process deduplicate {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.circular_contigs_deduplicate_memory}"
    def threads = "${params.circular_contigs_deduplicate_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
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
        awk -F "\t" '{OFS="\t"}  { if (\$1 != \$2) print \$0 }' ${blst_out_tsv} > ${id}/otherhit.tsv
        python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py ${circular_cut} > ${id}/circular_cut.all.length.txt
        ruby ${params.petagenomeDir}/scripts/Ruby/extract_contig_redundancy.3.rb -b ${id}/otherhit.tsv -l ${id}/circular_cut.all.length.txt -c ${params.circular_contigs_qc_thre_rd}  -d 6  -i ${id}/out_rd_info.txt  -o ${id}/out_ex_config.txt
        touch ${id}/circular.cut.fa ${id}/circular.extended.fa ${id}/circular.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_cut} | sed '/^\$/d'  > ${id}/circular.cut.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_ext} | sed '/^\$/d'  > ${id}/circular.extended.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular} | sed '/^\$/d' > ${id}/circular.fa
        """
}

// ==========================================
// 1. サブワークフロー（再利用可能な処理の本体）
// ==========================================

// 環状コンティグ検出・重複除去の一連処理
workflow CIRCULAR_CONTIGS_SUB {
    take:
    p
    contig

    main:
    if (params.use_pzlast) {
        // --- PZLAST モード ---
        if (params.test_pzlast_cfg1) {
            cfg1 = Channel.of(file(params.test_pzlast_cfg1))
        } else {
            cfg1 = Channel.of('')
        }

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
            tuple(p_val, ref_id, qry_id, tsv, cfg)
        }
        clsfy = classify(clsfy_p_in, contig)

        circular_cut = clsfy.map { id, cut, ext, circular, linear, selfhit_tsv ->
            tuple(id, cut)
        }

        // 3. PZLAST DB2 作成 ＋ アライメント 2
        if (params.test_pzlast_cfg2) {
            cfg2 = Channel.of(file(params.test_pzlast_cfg2))
        } else {
            cfg2 = Channel.of('')
        }

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
        ch_new = aln2.merge(clsfy).map {
            ref_id, qry_id, pzlst2_tsv, cfg,
            id, cut, ext, circular, linear, selfhit_tsv
            -> tuple(id, cut, ext, circular, pzlst2_tsv)
        }

        dedupl_in = p.combine(ch_new).map { p_val, id, cut, ext, circular, aln2_tsv ->
            tuple(p_val, id, cut, ext, circular, aln2_tsv)
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

        circular_cut = clsfy.map { id, cut, ext, circular, linear, selfhit_tsv ->
            tuple(id, cut)
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
        ch_new = aln2.merge(clsfy).map {
            ref_id, qry_id, blst2_tsv,
            id, cut, ext, circular, linear, selfhit_tsv
            -> tuple(id, cut, ext, circular, blst2_tsv)
        }

        dedupl_in = p.combine(ch_new).map { p_val, id, cut, ext, circular, aln2_tsv ->
            tuple(p_val, id, cut, ext, circular, aln2_tsv)
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
// 2. コマンドライン (-entry) 用エントリーポイント
// ==========================================

// A. メインの実行ワークフロー
workflow CIRCULAR_CONTIGS_ALL {
    p      = createNullParamsChannel()
    contig = createSeqsChannel(params.test_circular_contigs_contig)

    out_ch = CIRCULAR_CONTIGS_SUB(p, contig)

    // 各出力チャンネルの確認
    out_ch.aln1.view   { i -> "ALN1: $i" }
    out_ch.clsfy.view  { i -> "CLSFY: $i" }
    out_ch.aln2.view   { i -> "ALN2: $i" }
    out_ch.dedupl.view { i -> "DEDUPL: $i" }
}

// デフォルトエントリーポイント
workflow {
    CIRCULAR_CONTIGS_ALL()
}
