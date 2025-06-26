#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

workflow pool_contigs {
  take:
    reads
  main:
    merged = merge_contigs(contigs)
  emit:
    merged
}

workflow {
    contigs = channel.fromPath(params.test_pool_contigs_contigs, checkIfExists: true)
      .collect()
      .map{ it.sort() }
      .map{ it ->
         def key = it[0].simpleName.split('_')[0]
	 [key, it ] }
    contigs.view{ i -> "$i" }

//.map { contig_path -> tuple(contig_path.simpleName, contig_path) }
    //contigs = channel.fromFilePairs(params.test_pool_contigs_contigs, checkIfExists: true)
    //contigs.view{ i -> "$i" }
    //out = merge_contigs(contigs)
    //out.merged.view{ i -> "$i" }
}
