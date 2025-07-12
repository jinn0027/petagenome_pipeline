#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bbmap_bbmap_makerefdb_memory = params.memory
params.bbmap_bbmap_makerefdb_threads = params.threads

params.bbmap_bbmap_memory = params.memory
params.bbmap_bbmap_threads = params.threads

params.bbmap_ambiguous = "random"
params.bbmap_minid = 0.95
params.bbmap_pairlen = 1500

process bbmap_makerefdb {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/${task.process}", mode: 'copy', enabled: params.publish_output
    memory params.bbmap_bbmap_makerefdb_memory
    cpus params.bbmap_bbmap_makerefdb_threads

    input:
        tuple val(ref_id), path(ref, arity: '1')
    output:
        tuple val(ref_id), path("${ref_id}")
    script:
        """
        bbmap.sh \\
            -Xmx${params.bbmap_bbmap_makerefdb_memory}g \\
            threads=${params.bbmap_bbmap_threads} \\
            ref=${ref} \\
            path=${ref_id}
        """
}

process bbmap {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/${task.process}/${ref_id}", mode: 'copy', enabled: params.publish_output
    memory params.bbmap_bbmap_memory
    cpus params.bbmap_bbmap_threads

    input:
        tuple val(ref_id), path(ref_db, arity: '1'), val(pair_id), path(reads, arity: '2')
    output:
        tuple val(ref_id), val(pair_id), path("${pair_id}/out.sam", arity: '1')
    script:
        """
        mkdir -p ${pair_id}
        bbmap.sh \\
            -Xmx${params.bbmap_bbmap_memory}g \\
            threads=${params.bbmap_bbmap_threads} \\
            ambiguous=${params.bbmap_ambiguous} \\
            minid=${params.bbmap_minid} \\
            pairlen=${params.bbmap_pairlen} \\
            path=${ref_db} \\
            in=${reads[0]} \\
            in2=${reads[1]} \\
            out=${pair_id}/out.sam
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

