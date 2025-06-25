#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { spades_assembler } from "${params.petagenomeDir}/nf/spades"

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
//filter_contig.rename.py --min 1000 --rename --prefix {prefix}.{sample}.n. --table #{out}.1000.name.txt {dir_out}/scaffolds.fasta > {out}.1000.fa
//sample = a_samples[l]
//assembly_label = "SPAdes_meta"
//"#{dir_out_}/contig.#{sample}_#{assembly_label}"

process get_length {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    publishDir "${params.output}/get_length/${pair_id}", mode: 'copy'
    input:
        tuple val(pair_id), path(reads, arity: '1..*')

    output:
        tuple val(pair_id), path("out")
    script:
        """
        mkdir -p out
        reads_=( ${reads} )
        for i in \${reads_[@]}
        do
            zcat \$i | \
                awk '{if(\$1~/^\\+/||\$1~/^@/){print(\$1)}else{print(\$0)}}' | \
                python ${params.petagenomeDir}/scripts/Python/get_sequence_length.py -t fastq \
                >> out/\$(basename \$i).length.txt.gz
        done
        """
}

//filter_contig.rename.py --min 1000 --rename --prefix {prefix}.{sample}.n. --table #{out}.1000.name.txt {dir_out}/scaffolds.fasta > {out}.1000.fa
//filter_contig.rename.py --min 5000 {out}.1000.fa > {out}.5000.fa
//filter_contig.rename.py --min 10000 {out}.1000.fa > {out}.10000.fa


workflow {
   reads = channel.fromFilePairs(params.test_assembly_reads, checkIfExists: true)
   asm = spades_assembler(reads).map { pair_id, reads, unpaired -> tuple( pair_id, reads ) }
   asm.view { i -> "$i" }
   flt = filter_contig_rename(asm)
   flt.view { i -> "$i" }
   //fqc = fastqc(ec)
   //fqc.view{ i -> "$i" }
   //out = get_length(ec)
   //out.view{ i -> "$i" }
}
