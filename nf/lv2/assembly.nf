#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { printProcessProfile; createPairsChannel } from "${params.petagenomeDir}/nf/common/utils"
include { spades_assembler } from "${params.petagenomeDir}/nf/lv1/spades"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/lv1/blast"

params.assembly_filter_and_rename_memory = params.memory
params.assembly_filter_and_rename_threads = params.threads

params.assembly_get_length_memory = params.memory
params.assembly_get_length_threads = params.threads

params.assembly_get_stats_memory = params.memory
params.assembly_get_stats_threads = params.threads

process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "${params.apptainerRunOptions} --bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.assembly_filter_and_rename_memory} GB"
    cpus "${params.assembly_filter_and_rename_threads}"

    input:
        tuple val(id), path(read, arity: '1'), val(l_thre)
    output:
        tuple val(id), path("${id}/contig.${l_thre}.fa", arity: '0..*'), path("${id}/contig.name.txt")
    script:
        printProcessProfile(task)
        """
        mkdir -p ${id}
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min ${l_thre} --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.${l_thre}.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 1000 --rename --prefix n. --table ${id}/contig.name.txt ${read} > ${id}/contig.1000.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 5000 ${id}/contig.1000.fa > ${id}/contig.5000.fa
        #python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
        #     --min 10000 ${id}/contig.1000.fa > ${id}/contig.10000.fa
        """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.assembly_get_length_memory} GB"
    cpus "${params.assembly_get_length_threads}"

    input:
        tuple val(id), path(reads, arity: '0..*')
    output:
        tuple val(id), path("${id}/*.length.txt", arity: '0..*')
    script:
        printProcessProfile(task)
        """
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
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory "${params.assembly_get_stats_memory} GB"
    cpus "${params.assembly_get_stats_threads}"

    input:
        tuple val(id), path(lengths, arity: '1..*')
    output:
        tuple val(id), path("${id}/*.stats.txt")
    script:
        printProcessProfile(task)
        """
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

workflow assembly {
  take:
    reads
    l_thre
  main:
    asm = spades_assembler(reads)
    asm = asm.map { id, scaffolds, contigs ->
        tuple(id, 0 < scaffolds.size() ? scaffolds : contigs)
    }

    flt = filter_and_rename(asm.map { id, contigs -> tuple(id, contigs, l_thre) } )
    flt = flt.flatMap { id, contigs, name ->
        contigs.collect { c ->
	    if (c.size() != 0) {
              return [c.getBaseName(), c]
	    }
        }.findAll{ it != null }
    }

    len = get_length(flt)
    sts = get_stats(len)
    blstdb = blast_makerefdb(flt)

  emit:
    asm
    flt
    len
    sts
    blstdb
}

workflow {
    reads = createPairsChannel(params.test_assembly_reads)
    out = assembly(reads, params.test_assembly_l_thre)
    out.asm.view{ i -> "$i" }
    out.flt.view{ i -> "$i" }
    out.len.view{ i -> "$i" }
    out.sts.view{ i -> "$i" }
    out.blstdb.view{ i -> "$i" }
}
