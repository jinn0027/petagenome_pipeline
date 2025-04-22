#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.prodigal_procedure = "meta"
//params.prodigal_procedure = "single"
params.prodigal_outfmt = "gbk"
//params.prodigal_outfmt = "gff"
//params.prodigal_outfmt = "sco"

process prodigal {
    tag "${read_id}"
    container = "${params.petagenomeDir}/modules/prodigal/prodigal.sif"
    publishDir "${params.output}/prodigal/${read_id}", mode: 'copy'
    input:
        tuple val(read_id), path(read, arity: '1')
    output:
        tuple val(read_id), \
              path("out/${read_id}.faa", arity: '0..*'), \
              path("out/${read_id}.fna", arity: '0..*'), \
              path("out/${read_id}.${params.prodigal_outfmt}", arity: '0..*')
    script:
        """
        read_=${read}
        echo ${read} | grep -e ".gz\$" >& /dev/null && :
        if [ \$? -eq 0 ] ; then
            read_=\${read_%%.gz}
            unpigz -c ${read} > \${read_}
        fi
        mkdir -p out
        prodigal \\
            -p ${params.prodigal_procedure} \\
            -i \${read_} \\
            -f ${params.prodigal_outfmt} \\
            -a out/${read_id}.faa \\
            -d out/${read_id}.fna \\
            -o out/${read_id}.${params.prodigal_outfmt}
        """
}

workflow {
    read = channel.fromPath(params.test_prodigal_read, checkIfExists: true)
        .map { read_path -> tuple(read_path.simpleName, read_path) }
    out = prodigal(read)
    //out.view { i -> "${i}" }
}
