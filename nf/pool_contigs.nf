#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { cdhit_est } from "${params.petagenomeDir}/nf/cdhit"

//	#------------------------------
//	# concatenate assemblies
process merge_contigs {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(contigs, arity: '1..*')
    output:
        tuple val(id), path("out/merged_contig.fa")
    script:
        """
        mkdir -p out
        touch out/merged_contig.fa
        contigs_=( ${contigs} )
        for contig in \${contigs_[@]}
        do
            contig_=\${contig}
            echo \${contig} | grep -e ".gz\$" >& /dev/null && :
            if [ \$? -eq 0 ] ; then
                contig_=\${contig_%%.gz}
                unpigz -c \${contig} > \${contig_}
            fi
            cat \${contig_} >> out/merged_contig.fa
        done
        """
}

//	# rename and filter (L>=#{l_thre}) contigs
//	# this contig name format is nessecary for later procedures.
//	# the corresponding name table will be output in #{out}.#{l_thre}.name.txt.
//	f_qsub.puts "python ${SCRIPT_RENAME_} -min #{l_thre} --rename --prefix #{prefix}.#{project}.#{merge_label}.n. --table #{out}.#{l_thre}.name.txt #{rep_fa_} > #{out}.#{l_thre}.fa"


//	# summaize contig name ( contig in each sample / representative contig / renamed representative contig )
//	f_qsub.puts "ruby ${SCRIPT_CLSTR_} -i #{rep_fa_}.clstr --include_rep > #{rep_fa_}.name"
//	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{rep_fa_}.name 2 #{out}.#{l_thre}.name.txt 2 | awk -F '\\t' '{OFS=\"\\t\"} {print $2,$1,$3}' > #{out_table_name_}"
process parse_cdhit_cruster {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(merged, arity: '1')
    output:
        tuple val(id), path("out/${id}.name")
    script:
        """
        mkdir -p out
        echo ruby ${params.petagenomeDir}/scripts/Ruby/parse.cdhit_clstr.rb -i ${merged} --include_rep > out/${id}.name
        """
}

workflow pool_contigs {
  take:
    contigs
  main:
    // concatenate assemblies
    merged = merge_contigs(contigs)

    // removing redundancy by CD-HIT-EST
    cdhit = cdhit_est(merged)

    //
    cdhit_ = cdhit.map{ id, fasta, clstr -> tuple(id, fasta) }

    // rename and filter (L>=#{l_thre}) contigs
    // this contig name format is nessecary for later procedures.
    // the corresponding name table will be output in #{out}.#{l_thre}.name.txt.

    //name = parse_cdhit_cruster(cdhit_)
  emit:
    merged
    cdhit
    //name
}

workflow {
    contigs = channel.fromPath(params.test_pool_contigs_contigs, checkIfExists: true)
      .collect()
      .map{ it.sort() }
      .map{ it ->
         def key = it[0].simpleName.split('_')[0]
	 [key, it ] }
    contigs.view{ i -> "$i" }
    out = pool_contigs(contigs)
    //out.merged.view{ i -> "$i" }
    out.cdhit.view{ i -> "$i" }
    //out.name.view{ i -> "$i" }
}
