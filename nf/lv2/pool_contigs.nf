#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { clusterOptions; processProfile; createSeqsChannel } from "${params.petagenomeDir}/nf/common/utils"
include { cdhit_est } from "${params.petagenomeDir}/nf/lv1/cdhit"
include { mmseqs2_makerefdb; mmseqs2_cluster } from "${params.petagenomeDir}/nf/lv1/mmseqs2"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"

params.pool_contigs_mergs_contigs_memory = params.memory
params.pool_contigs_mergs_contigs_threads = params.threads

params.pool_contigs_filter_and_rename_memory = params.memory
params.pool_contigs_filter_and_rename_threads = params.threads

params.pool_contigs_summarize_name_memory = params.memory
params.pool_contigs_summarize_name_threads = params.threads

params.pool_contigs_get_length_memory = params.memory
params.pool_contigs_get_length_threads = params.threads

params.pool_contigs_get_stats_memory = params.memory
params.pool_contigs_get_stats_threads = params.threads

//params.pool_contigs_clustering_process = "cdhit"
params.pool_contigs_clustering_process = "mmseqs2"

process merge_contigs {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.pool_contigs_mergs_contigs_memory} GB"
    threads = "${params.pool_contigs_mergs_contigs_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(contigs, arity: '1..*', stageAs: "?/*")
    output:
        tuple val(id), path("${id}/merged_contig.fa"), path("${id}/contig_list.txt")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        touch ${id}/merged_contig.fa
        touch ${id}/contig_list.txt
        contigs_=( ${contigs} )
        for contig in \${contigs_[@]}
        do
            contig_=\${contig}
            echo \${contig} | grep -e ".gz\$" >& /dev/null && :
            if [ \$? -eq 0 ] ; then
                contig_=\${contig_%%.gz}
                unpigz -c \${contig} > \${contig_}
            fi
            cat \${contig_} >> ${id}/merged_contig.fa
            echo \${contig_} >> ${id}/contig_list.txt
        done
        """
}

process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.pool_contigs_filter_and_rename_memory} GB"
    threads = "${params.pool_contigs_filter_and_rename_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(read, arity: '1'), val(l_thre)
    output:
        tuple val(id), path("${id}/contig.${l_thre}.fa"), path("${id}/contig.name.txt")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min ${l_thre} --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.${l_thre}.fa
        """
}

process summarize_name {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.pool_contigs_summarize_name_memory} GB"
    threads = "${params.pool_contigs_summarize_name_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(name, arity: '1')
        tuple val(id), path(clstr, arity: '1')
    output:
        tuple val(id), path("${id}/${id}.name.txt"), path("${id}/*")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        if [ "${params.pool_contigs_clustering_process}" = "mmseqs2" ] ; then
            cp -f ${clstr} ${id}/${id}.name_
        else
            ruby ${params.petagenomeDir}/scripts/Ruby/parse.cdhit_clstr.rb -i ${clstr} --include_rep > ${id}/${id}.name_
        fi
        ruby ${params.petagenomeDir}/scripts/Ruby/join_with_tab.rb ${id}/${id}.name_ 2 ${name} 2 | \
             awk -F '\\t' '{OFS=\"\\t\"} {print \$2,\$1,\$3}' > ${id}/${id}.name.txt
        awk '{print(\$1)}' ${id}/${id}.name.txt | sort | uniq > ${id}/${id}.samples.txt
        while read sample
        do
            awk -F \"\\t\" -v sample=\${sample} \
                '{OFS=\"\\t\"} { if (\$1 ~ sample) print \$0 }' \
                ${id}/${id}.name.txt \
                > ${id}/${id}.\${sample}.name.txt
        done<${id}/${id}.samples.txt
        rm -f ${id}/${id}.name_
        """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.pool_contigs_get_length_memory} GB"
    threads = "${params.pool_contigs_get_length_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(reads, arity: '1..*')
    output:
        tuple val(id), path("${id}/*.length.txt")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' \${i} | \
            python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fasta \
            > ${id}/\$(basename \$i).length.txt
        done
        """
}

process get_stats {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.pool_contigs_get_stats_memory} GB"
    threads = "${params.pool_contigs_get_stats_threads}"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, threads, label)}"
    input:
        tuple val(id), path(lengths, arity: '1..*')
    output:
        tuple val(id), path("${id}/*.stats.txt")
    script:
        """
        echo "${processProfile(task)}"
        mkdir -p ${id}
        lengths_=( ${lengths} )
        for i in \${lengths_[@]}
        do
          outname=\$(basename \${i} | sed "s#length#stats#")
          n=\$(cat \${i} | wc -l)
          if [ \${n} -gt 0 ] ; then
              #R --vanilla --slave --args \${i} 2 <  ${params.petagenomeDir}/scripts/R/stats.assembly.R > ${id}/\${outname}
              Rscript ${params.petagenomeDir}/scripts/R/stats.assembly.R \${i} 2 > ${id}/\${outname}
          else
              touch ${id}/\${outname}
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

    // removing redundancy
    if ( params.pool_contigs_clustering_process == "mmseqs2" ) {
        merged_db = mmseqs2_makerefdb( merged.map{ id, fasta, list -> tuple(id, fasta ) } )
        clust = mmseqs2_cluster(merged_db)
    } else {
        clust = cdhit_est( merged.map{ id, fasta, list -> tuple(id, fasta ) } )
    }

    // rename and filter (L>=#{l_thre}) contigs
    flt = filter_and_rename( clust.map{ id, fasta, clstr -> tuple(id, fasta, l_thre) } )

    // summarize contig name ( contig in each sample / representative contig / renamed representative contig )
    name = summarize_name( flt.map{ id, fasta, name -> tuple(id, name) }, clust.map{ id, fasta, clstr -> tuple(id, clstr) } )

    // get length of contigs
    len = get_length( flt.map{ id, fasta, name -> tuple(id, fasta) } )

    // stats of assemblies
    sts = get_stats(len)

    // blastdb
    blstdb = blast_makerefdb( flt.map{ id, contig, name -> tuple(id, contig) } )

  emit:
    merged
    clust
    flt
    name
    len
    sts
    blstdb
}

workflow {
    def contigs = createSeqsChannel(params.test_pool_contigs_contigs)
    contigs.view{ i -> "$i" }
    out = pool_contigs(contigs, params.test_pool_contigs_l_thre)
    //out.merged.view{ i -> "$i" }
    //out.clust.view{ i -> "$i" }
    //out.flt.view{ i -> "$i" }
    //out.name.view{ i -> "$i" }
    //out.len.view{ i -> "$i" }
    //out.sts.view{ i -> "$i" }
    out.blstdb.view{ i -> "$i" }
}
