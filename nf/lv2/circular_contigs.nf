#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.circular_contigs_select_selfhit_memory = params.memory
params.circular_contigs_select_selfhit_threads = params.threads

// E-value cutoff for circular formation
params.circular_contigs_e_thre = "1e-10"
// % identity for circular formation
params.circular_contigs_pi_thre = "100"
// minimum alignment length for circular formation
params.circular_contigs_al_thre = "50"
// % identity for removal of redundancy
params.circular_contigs_pi_thre_rd = "95"
// % query coverage for removal of redundancy
params.circular_contigs_qc_thre_rd = "95"
// minimum length for linear contigs
params.circular_contigs_len_l = "5000"
// minimum length for circular contigs
params.circular_contigs_len_c = "1500"

params.test_circular_contigs_l_thre = 1000
//params.test_circular_contigs_l_thre = 5000

include { createNullParamsChannel; getParam; clusterOptions; processProfile } \
    from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb; blastn } from "${params.petagenomeDir}/nf/lv1/blast"

process select_selfhit {
    tag "${ref_id}_@_${qry_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
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
    output:
        tuple val(ref_id), val(qry_id), path("${qry_id}/selfhit.tsv", arity: '1')
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${qry_id}
        awk -F "\t" '{OFS="\t"}  { if (\$1 == \$2) print \$0 }' ${in_tsv} > ${qry_id}/selfhit.tsv
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

    //${BLASTN_} -task #{task} -num_threads #{n_threads} -query #{in_contig_} \
    //           -db #{in_contig} -evalue #{e_thre} -perc_identity #{pi_self} -outfmt 6 -num_alignments 5 -out #{out_blast_all_}
    blstn = blastn(p, blstin)
    selfhit = select_selfhit(p, blstn)
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
