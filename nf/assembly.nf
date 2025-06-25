#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { spades_assembler } from "${params.petagenomeDir}/nf/spades"
include { blast_makerefdb } from "${params.petagenomeDir}/nf/blast"

process filter_contig_rename {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy'
    input:
        tuple val(id), path(read, arity: '1')
    output:
        tuple val(id), path("out/*.fa", arity:"3")
    script:
        def assembly_label = "SPAdes_meta"
        def prefix
        if (params.assembly_type == "virome") {
            prefix = "v"
        } else if (params.assembly_type == "bacteriome") {
            previx = "b"
        } else {
            throw new Exception("Error: illegal assemlby_type ${params.assembly_type}")
        }
        """
        mkdir -p out
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 1000 --rename --prefix ${prefix}.${id}.n. --table out/contig.${id}_${assembly_label}.name.txt ${read} > out/contig.${id}_${assembly_label}.1000.fa
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 5000 out/contig.${id}_${assembly_label}.1000.fa > out/contig.${id}_${assembly_label}.5000.fa
        python ${params.petagenomeDir}/scripts/Python/filter_contig.rename.py \
             --min 10000 out/contig.${id}_${assembly_label}.1000.fa > out/contig.${id}_${assembly_label}.10000.fa
        """
}

process get_length {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy'
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
            cat \$i | \
                awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' | \
                python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fasta \
                > out/\$(basename \$i).length.txt
        done
        """
}

process stats_assembly {
    tag "${id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${id}", mode: 'copy'
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
    //asm.view { i -> "$i" }
    flt = filter_contig_rename(asm)
    //flt.view { i -> "$i" }
    //log.info "${flt.dump()}"

    len = get_length(flt)
    //len.view{ i -> "$i" }
    sts = stats_assembly(len)
    //sts.view{ i -> "$i" }

    refs = flt.flatMap { key, files ->
        files.collect{ file_path ->
	    if (file_path.size() != 0) {
              return [file_path.getBaseName(), file_path]
	    }
        }.findAll{ it != null }
    }
    //refs.view( i -> "$i" )

    blstdb = blast_makerefdb(refs)
    //blstdb.view{ i -> "$i" }

  emit:
    asm
    flt
    len
    sts
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
