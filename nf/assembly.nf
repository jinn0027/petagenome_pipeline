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
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '1..*')
    output:
        tuple val(pair_id), path("out/*.length.txt")
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
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/${task.process}/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(lengths, arity: '1..*')
    output:
        tuple val(pair_id), path("out/*.stats.txt")
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

workflow {
   reads = channel.fromFilePairs(params.test_assembly_reads, checkIfExists: true)
   asm = spades_assembler(reads).map { pair_id, reads, unpaired -> tuple( pair_id, reads ) }
   asm.view { i -> "$i" }
   flt = filter_contig_rename(asm)
   //flt.view { i -> "$i" }
   len = get_length(flt)
   len.view{ i -> "$i" }
   sts = stats_assembly(len)
   sts.view{ i -> "$i" }
   params.blast_dbtype = "nucl"
   //fqc = fastqc(ec)
   //fqc.view{ i -> "$i" }
}
