#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { printProcessProfile } from "${params.petagenomeDir}/nf/common/utils"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"

params.circular_contigs_explore_circular_contigs_memory = params.memory
params.circular_contigs_explore_circular_contigs_threads = params.threads

//# E-value cutoff for circular formation
//e_thre = "1e-10"
//# % identity for circular formation
//pi_thre = "100"
//# minimum alignment length for circular formation
//al_thre = "50"
//
//# % identity for removal of redundancy
//pi_thre_rd = "95"
//# % query coverage for removal of redundancy
//qc_thre_rd = "95"
//
//# minimum length for linear contigs
//len_l = "5000"
//# minimum length for circular contigs
//len_c = "1500"

//	# main script for exploring circular contigs
//	ruby ${SCRIPT_SELF_}  -d #{dir_out_}  -i #{query_fa_}  -n #{n_threads}  --len_l #{len_l} --len_c #{len_c} --pi_rd #{pi_thre_rd} --qc_rd #{qc_thre_rd} "

process explore_circular_contigs {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    memory "${params.circular_contigs_explore_circular_contigs_memory} GB"
    cpus "${params.circular_contigs_explore_circular_contigs_threads}"

    def 

    input:
        tuple val(id), path(contig, arity: '1')
    output:
        tuple val(id), path("out/*")
    script:
        printProcessProfile(task)
        """
        mkdir -p out
        """
}


//	# copy results to adjust to the existing scripts
//	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.fa  #{out}.fa"
//	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.length.txt  #{out}.length.txt"
//	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.extended.fa  #{out}.extended.fa"
//	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.extended.length.txt  #{out}.extended.length.txt"


//# blastdb
//	f_qsub.puts "${MAKEBLASTDB_} -in #{out}.fa -out #{out} -dbtype nucl -parse_seqids"

//	# name information
//	# updated (final) contigs
//	f_qsub.puts "cut -f 1 #{out}.length.txt > #{out}.txt"
//	f_qsub.puts "cat #{out}.txt | perl -pe \"s/l(\\.\\d+)$/n\\1/\" | perl -pe \"s/c(\\.\\d+)$/n\\1/\" | paste - #{out}.txt > #{dir_out_}/name.updated.#{query_label}.#{len_l}.c_#{len_c}.txt"
//	# redundant (excluded) contigs
//	f_qsub.puts "cut -f 1,6 #{dir_out_}/contig.redundant.c_#{len_c}.length.txt | perl -pe \"s/l(\\.\\d+)\\t/n\\1\\t/\" | perl -pe \"s/c(\\.\\d+)\\t/n\\1\\t/\" > #{dir_out_}/name.redundant.#{query_label}.#{len_l}.c_#{len_c}.txt"
//	# all contigs
//	f_qsub.puts "cat #{dir_out_}/name.updated.#{query_label}.#{len_l}.c_#{len_c}.txt #{dir_out_}/name.redundant.#{query_label}.#{len_l}.c_#{len_c}.txt > #{dir_out_}/name.all.#{query_label}.#{len_l}.c_#{len_c}.txt"


workflow circular_contigs {
  take:
    contig
  main:
    
  emit:
    contig
}

workflow {
    def l_thre = "1000" // virome
    //def l_thre = "5000" // bacteriome

    contig = channel.fromPath(params.test_circular_contigs_contig, checkIfExists: true)
      .map{ path -> tuple(path.simpleName, path) }
    contig.view { i -> "$i" }

    out = circular_contigs(contig)
    out.contig.view { i -> "$i" }
}
