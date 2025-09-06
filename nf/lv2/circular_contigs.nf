#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile } from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb; blastn } from "${params.petagenomeDir}/nf/lv1/blast"

params.circular_contigs_explore_circular_contigs_memory = params.memory
params.circular_contigs_explore_circular_contigs_threads = params.threads

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

workflow circular_contigs {
  take:
    contig
    l_thre
  main:
    ref = contig
    qry = contig
    ref_db = blast_makerefdb(ref)
    in = ref_db.combine(qry)
    out = blastn(in)
  emit:
    out
}

workflow {
    contig = channel.fromPath(params.test_circular_contigs_contig, checkIfExists: true)
      .map{ path -> tuple(path.simpleName, path) }
    contig.view { i -> "$i" }
    out = circular_contigs(contig, params.test_circular_contigs_l_thre)
    out.out.view { i -> "$i" }
}
