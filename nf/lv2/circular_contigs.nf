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
params.test_circular_contigs_l_thre = 1000
//params.test_circular_contigs_l_thre = 5000

include { createNullParamsChannel; getParam; clusterOptions; processProfile } \
    from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb as blast_makerefdb1; blastn as blastn1} from "${params.petagenomeDir}/nf/lv1/blast"
include { blast_makerefdb as blast_makerefdb2; blastn as blastn2} from "${params.petagenomeDir}/nf/lv1/blast"

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
        val(p)
        tuple val(ref_id), val(qry_id), path(in_tsv, arity: '1')
        tuple val(qry_id), path(in_qry, arity: '1')
    output:
        tuple val(ref_id),
              path("${qry_id}/circular.cut.fa"),
              path("${qry_id}/circular.extended.fa"),
              path("${qry_id}/circular.fa"),
              path("${qry_id}/linear.fa"),
              path("${qry_id}/selfhit.tsv", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        awk -F "\t" '{OFS="\t"}  { if (\$1 == \$2) print \$0 }' ${in_tsv} > ${qry_id}/selfhit.tsv
        touch ${qry_id}/circular.cut.fa ${qry_id}/circular.extended.fa ${qry_id}/circular.fa ${qry_id}/linear.fa 
        seqkit fx2tab -j ${params.circular_contigs_classify_threads} -n -i -l ${in_qry} | while read -r id len; do
            pos_end=\$(awk -v id=\${id} -v len=\${len} -v al_self=${getParam(p, 'circular_contigs_al_self')} \\
                '{if (\$1==id && \$4!=len && \$4>=al_self && \$9==1) print(\$7-1)}' \\
                ${qry_id}/selfhit.tsv | sort -r | head -n 1)
            seqkit grep -np \${id} ${in_qry} >> _fa
            if [ "\${pos_end}" != "" ] ; then
                cat _fa | seqkit subseq -r 1:\${pos_end} >> ${qry_id}/circular.cut.fa
                cp ${qry_id}/circular.cut.fa ${qry_id}/circular.extended.fa
                cat _fa | seqkit seq -s >> ${qry_id}/circular.extended.fa
                cat _fa >> ${qry_id}/circular.fa
            else
                cat _fa >> ${qry_id}/linear.fa
            fi
            rm -f _fa
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
        val(p)
        tuple val(id), path(circular_cut), path(circular_ext), path(circular), path(blst_out_tsv, arity: '1')
    output:
        tuple val(id),
              path("${id}/circular.cut.fa"),
              path("${id}/circular.extended.fa"),
              path("${id}/circular.fa"),
              path("${id}/otherhit.tsv", arity: '1'),
              path("${id}/*.txt", arity: '3')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        awk -F "\t" '{OFS="\t"}  { if (\$1 != \$2) print \$0 }' ${blst_out_tsv} > ${id}/otherhit.tsv
        python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py ${circular_cut} > ${id}/circular_cut.all.length.txt
        ruby ${params.petagenomeDir}/scripts/Ruby/extract_contig_redundancy.3.rb -b ${id}/otherhit.tsv -l ${id}/circular_cut.all.length.txt -c ${params.circular_contigs_qc_thre_rd}  -d 6  -i ${id}/out_rd_info.txt  -o ${id}/out_ex_config.txt
        touch ${id}/circular.cut.fa ${id}/circular.extended.fa ${id}/circular.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_cut} > ${id}/circular.cut.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular_ext} > ${id}/circular.extended.fa
        python ${params.petagenomeDir}/scripts/Python/filter_fasta_by_id.py -f ${id}/out_ex_config.txt ${circular} > ${id}/circular.fa
        """
}

workflow circular_contigs {
  take:
    p
    contig
    l_thre
  main:
    blstdb1 = blast_makerefdb1(p, contig)
    blstin1 = blstdb1.combine(contig)
    p_blastn1 = Channel.of(['blast_task':'megablast',
                            'blast_perc_identity':params.circular_contigs_pi_self,
                            'blast_evalue':params.circular_contigs_e_thre,
                            'blast_outfmt':6,
                            'blast_num_alignments':params.circular_contigs_blast1_num_alignments
                            ])
    blstn1 = blastn1(p_blastn1, blstin1)
    clsfy = classify(p, blstn1, contig)

    circular_cut = clsfy.map { id, circular_cut, circular_extended, circular, linear, selfhit_tsv ->
        [ id, circular_cut ]
    }

    blstdb2 = blast_makerefdb2(p, circular_cut)
    p_blastn2 = Channel.of(['blast_task':'megablast',
                            'blast_perc_identity':params.circular_contigs_pi_self,
                            'blast_evalue':params.circular_contigs_e_thre,
                            'blast_outfmt':6,
                            'blast_num_alignments':params.circular_contigs_blast2_num_alignments
                            ])
    blstin2 = blstdb1.combine(circular_cut)
    blstn2 = blastn2(p_blastn2, blstin2)

    blstn2.view{ i-> "$i" }
    clsfy.view{ i-> "$i" }
    ch_new = blstn2.merge(clsfy).map {
        ref_id, qry_id, blst2_tsv,
        id, cut, ext, circular, linear, selfhit_tsv
        -> [id, cut, ext, circular, blst2_tsv]
    }
    ch_new.view{ i -> "${i}" }
    dedupl = deduplicate(p, ch_new)
  emit:
    blstn1
    blstn2
    clsfy
    dedupl
}

workflow {
    p = createNullParamsChannel()
    contig = channel.fromPath(params.test_circular_contigs_contig, checkIfExists: true)
      .map{ path -> tuple(path.simpleName, path) }
    contig.view { i -> "$i" }
    out = circular_contigs(p, contig, params.test_circular_contigs_l_thre)
    //out.blstn1.view { i -> "$i" }
    //out.clsfy.view { i -> "$i" }
}
