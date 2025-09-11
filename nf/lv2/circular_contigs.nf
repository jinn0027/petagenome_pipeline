#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.circular_contigs_select_selfhit_memory = params.memory
params.circular_contigs_select_selfhit_threads = params.threads

params.circular_contigs_e_thre = "1e-10"        // E-value cutoff for circular formation
params.circular_contigs_pi_thre = "100"         // identity for circular formation
params.circular_contigs_al_thre = "50"          // minimum alignment length for circular formation
params.circular_contigs_pi_thre_rd = "95"       // identity for removal of redundancy
params.circular_contigs_qc_thre_rd = "95"       // query coverage for removal of redundancy
params.circular_contigs_len_l = "5000"          // minimum length for linear contigs
params.circular_contigs_len_c = "1500"          // minimum length for circular contigs
params.circular_contigs_pi_self = 100           // identity for circular formation
params.circular_contigs_al_self = 50            // alignment length for circular formation
params.circular_contigs_blast_num_alignments=5
params.test_circular_contigs_l_thre = 1000
//params.test_circular_contigs_l_thre = 5000

include { createNullParamsChannel; getParam; clusterOptions; processProfile } \
    from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb; blastn } from "${params.petagenomeDir}/nf/lv1/blast"

process select_selfhit {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/seqkit/seqkit.sif"
    //containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    def gb = "${params.circular_contigs_select_selfhit_memory}"
    def threads = "${params.circular_contigs_select_selfhit_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    input:
        val(p)
        tuple val(ref_id), val(qry_id), path(in_tsv, arity: '1')
        tuple val(qry_id), path(in_qry, arity: '1')
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/selfhit.tsv", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        awk -F "\t" '{OFS="\t"}  { if (\$1 == \$2) print \$0 }' ${in_tsv} > ${qry_id}/selfhit.tsv
        seqkit stats -j 16 ${in_qry} >> ${qry_id}/selfhit.tsv
        """
}

workflow circular_contigs {
  take:
    p
    contig
    l_thre
  main:
    ref = contig
    qry = contig
    blstdb = blast_makerefdb(p, ref)
    blstin = blstdb.combine(qry)
    p_blastn = Channel.of(['blast_task':'megablast',
                           'blast_perc_identity':params.circular_contigs_pi_self,
                           'blast_evalue':params.circular_contigs_e_thre,
                           'blast_outfmt':6,
                           'blast_num_alignments':params.circular_contigs_blast_num_alignments
                           ])
    blstn = blastn(p_blastn, blstin)
    selfhit = select_selfhit(p, blstn, qry)
  emit:
    blstn
    selfhit
}

workflow {
    p = createNullParamsChannel()
    contig = channel.fromPath(params.test_circular_contigs_contig, checkIfExists: true)
      .map{ path -> tuple(path.simpleName, path) }
    contig.view { i -> "$i" }
    out = circular_contigs(p, contig, params.test_circular_contigs_l_thre)
    out.blstn.view { i -> "$i" }
    out.selfhit.view { i -> "$i" }
}
