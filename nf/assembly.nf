#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { spades_assembler } from "${params.petagenomeDir}/nf/spades"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/blast"

process filter_and_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(read, arity: '1')
    output:
        tuple val(id), path("out/contig.*.fa", arity:"3"), path("out/contig.name.txt")
    script:
        """
        mkdir -p out
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 1000 --rename --prefix n. --table out/contig.name.txt ${read} > out/contig.1000.fa
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 5000 out/contig.1000.fa > out/contig.5000.fa
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 10000 out/contig.1000.fa > out/contig.10000.fa
        """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(id), path(reads, arity: '3')
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

process get_stats {
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

workflow assembly {
  take:
    reads
  main:
    asm = spades_assembler(reads).map { id, paired, unpaired -> tuple( id, paired ) }
    flt = filter_and_rename(asm)
    len = get_length(flt.map{id, contigs, name -> tuple(id, contigs)})
    sts = get_stats(len)
    ctg = flt.flatMap { id, contigs, name ->
        contigs.collect{ contig_path ->
	    if (contig_path.size() != 0) {
              return [contig_path.getBaseName(), contig_path]
	    }
        }.findAll{ it != null }
    }
    blstdb = blast_makerefdb(ctg)
  emit:
    asm
    flt
    len
    sts
    ctg
    blstdb
}

workflow {
    reads = channel.fromFilePairs(params.test_assembly_reads, checkIfExists: true)
    out = assembly(reads)
    out.asm.view{ i -> "$i" }
    out.flt.view{ i -> "$i" }
    out.len.view{ i -> "$i" }
    out.sts.view{ i -> "$i" }
    out.blstdb.view{ i -> "$i" }
}
