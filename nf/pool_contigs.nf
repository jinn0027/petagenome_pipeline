#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { cdhit_est } from "${params.petagenomeDir}/nf/cdhit"

// concatenate assemblies
process merge_contigs {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(contigs, arity: '1..*')
    output:
        tuple val(id), path("out/merged_contig.fa"), path("out/contig_list.txt")
    script:
        """
        mkdir -p out
        touch out/merged_contig.fa
        touch out/contig_list.txt
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
            echo \${contig_} >> out/contig_list.txt
        done
        """
}

// rename and filter (L>=#{l_thre}) contigs
process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(read, arity: '1'), val(l_thre)
    output:
        tuple val(id), path("out/contig.${l_thre}.fa"), path("out/contig.name.txt")
    script:
        """
        mkdir -p out
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min ${l_thre} --rename --prefix n. --table out/contig.name.txt ${read} > out/contig.${l_thre}.fa
        """
}

// summarize contig name ( contig in each sample / representative contig / renamed representative contig )
process summarize_name {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(name, arity: '1')
        tuple val(id), path(clstr, arity: '1')
    output:
        tuple val(id), path("out/${id}.name.txt"), path("out/*")
    script:
        """
        mkdir -p out
        ruby ${params.petagenomeDir}/scripts/Ruby/parse.cdhit_clstr.rb -i ${clstr} --include_rep > out/${id}.name_
        ruby ${params.petagenomeDir}/scripts/Ruby/join_with_tab.rb out/${id}.name_ 2 ${name} 2 | \
             awk -F '\\t' '{OFS=\"\\t\"} {print \$2,\$1,\$3}' > out/${id}.name.txt
        awk '{print(\$1)}' out/${id}.name.txt | sort | uniq > out/${id}.samples.txt
        while read sample
        do
            awk -F \"\\t\" -v sample=\${sample} '{OFS=\"\\t\"} { if (\$1 ~ sample) print \$0 }' out/${id}.name.txt > out/${id}.\${sample}.name.txt
        done<out/${id}.samples.txt
        """
}

//	f_qsub.puts "R --vanilla --slave --args #{out}.#{l_thre}.length.txt 2 < ${SCRIPT_STATS_} > #{out}.#{l_thre}.stats.txt"
//	# blastdb
//	f_qsub.puts "${MAKEBLASTDB_} -in #{out}.#{l_thre}.fa -out #{out}.#{l_thre} -dbtype nucl -parse_seqids"

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(reads, arity: '1..*')
    output:
        tuple val(id), path("out/*.length.txt")
    script:
        """
        mkdir -p out
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' \${i} | \
            python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fasta \
            > out/\$(basename \$i).length.txt
        done
        """
}

process stats_assembly {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(lengths, arity: '1..*')
    output:
        tuple val(id), path("out/*.stats.txt")
    script:
        """
        mkdir -p out
        lengths_=( ${lengths} )
        for i in \${lengths_[@]}
        do
          outname=\$(basename \${i} | sed "s#length#stats#")
          n=\$(cat \${i} | wc -l)
          if [ \${n} -gt 0 ] ; then
              #R --vanilla --slave --args \${i} 2 <  ${params.petagenomeDir}/scripts/R/stats.assembly.R > out/\${outname}
              Rscript ${params.petagenomeDir}/scripts/R/stats.assembly.R \${i} 2 > out/\${outname}
          else
              touch out/\${outname}
          fi
        done
        """
}

workflow pool_contigs {
  take:
    contigs
    l_thre
  main:
    // concatenate assemblies
    merged = merge_contigs(contigs)

    // removing redundancy by CD-HIT-EST
    cdhit = cdhit_est( merged.map{ id, fasta, list -> tuple(id, fasta ) } )
    cdhit_ = cdhit.map{ id, fasta, clstr -> tuple(id, fasta, l_thre) }

    cdhit_.view{ i -> "$i" } 

    // rename and filter (L>=#{l_thre}) contigs
    flt = filter_and_rename(cdhit_)

    // summarize contig name ( contig in each sample / representative contig / renamed representative contig )
    name = summarize_name(
        flt.map{ id, fasta, name -> tuple(id, name) },
	cdhit.map{ id, fasta, clstr -> tuple(id, clstr) }
    )
    name.view{ i -> "$i" } 

    len = get_length( flt.map{ id, fasta, name -> tuple(id, fasta) } )
    sts = stats_assembly(len)

  emit:
    merged
    cdhit
    flt
    name
    len
    sts
}

workflow {
    def l_thre = "1000"

    contigs = channel.fromPath(params.test_pool_contigs_contigs, checkIfExists: true)
      .collect()
      .map{ it.sort() }
      .map{ it ->
         def key = it[0].simpleName.split('_')[0]
	 [key, it ] }
    contigs.view{ i -> "$i" }
    out = pool_contigs(contigs, l_thre)
    out.merged.view{ i -> "$i" }
    //out.cdhit.view{ i -> "$i" }
    //out.flt.view{ i -> "$i" }
    //out.name.view{ i -> "$i" }
}
