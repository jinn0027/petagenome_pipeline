#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bbmap_ambiguous = "random"
params.bbmap_minid = 0.95
params.bbmap_pairlen = 1500

process bbmap_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("ref_db")
    script:
        """
        bbmap.sh \\
            -Xmx${params.memory}g \\
            threads=${params.threads} \\
            ref=${ref} \\
            path=ref_db
        """
}

process bbmap {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/${task.process}/${ref_id}/${pair_id}", mode: 'copy', enabled: params.publish_output
    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(ref_id), val(pair_id), path("out.sam", arity: '1')
    script:
        """
        mkdir -p out
        bbmap.sh \\
            -Xmx${params.memory}g \\
            threads=${params.threads} \\
            ambiguous=${params.bbmap_ambiguous} \\
            minid=${params.bbmap_minid} \\
            pairlen=${params.bbmap_pairlen} \\
            path=${ref_db} \\
            in=${reads[0]} \\
            in2=${reads[1]} \\
            out=out.sam
        """
}

workflow {
    ref = channel.fromPath(params.test_bbmap_ref, checkIfExists: true)
        .map { ref_path -> tuple(ref_path.simpleName, ref_path) }
    reads = channel.fromFilePairs(params.test_bbmap_reads, checkIfExists: true)

    //ref.view { i -> "$i" }
    //reads.view { i -> "$i" }

    ref_db = bbmap_makerefdb(ref)
    //ref_db.view { i -> "$i" }
    in = ref_db.combine(reads)
    //in.view { i -> "$i" }
    out = bbmap(in)
    out.view { i -> "$i" }
}

