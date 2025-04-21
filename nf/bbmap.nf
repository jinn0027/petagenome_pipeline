#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bbmap_ambiguous = "random"
params.bbmap_minid = 0.95
params.bbmap_pairlen = 1500

params.test_bbmap_separate = true

process bbmap_make_index {
    tag "${ref_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/bbmap/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref)

    output:
        tuple val(ref_id), path("index")
    script:
    """
    bbmap.sh \\
        -Xmx${params.memory}g \\
        threads=${params.threads} \\
        ref=${ref} \\
        path=index
    """
}

process bbmap_align {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/bbmap/${ref_id}/${pair_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(index), val(pair_id), path(reads)

    output:
        tuple val(ref_id), val(pair_id), path("*.sam")
    script:
    """
    bbmap.sh \\
        -Xmx${params.memory}g \\
        threads=${params.threads} \\
        ambiguous=${params.bbmap_ambiguous} \\
        minid=${params.bbmap_minid} \\
        pairlen=${params.bbmap_pairlen} \\
        path=${index} \\
        in=${reads[0]} \\
        in2=${reads[1]} \\
        out=${ref_id}_@_${pair_id}_bbmap_out.sam
    """
}

process bbmap {
    tag "${ref_id}_@_${pair_id}"
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/bbmap/${ref_id}/${pair_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref), val(pair_id), path(reads)

    output:
        tuple val(ref_id), val(pair_id), path("*.sam")
    script:
    """
    bbmap.sh \\
        -Xmx${params.memory}g \\
        threads=${params.threads} \\
        ref=${ref} \\
        path=index
    bbmap.sh \\
        -Xmx${params.memory}g \\
        threads=${params.threads} \\
        ambiguous=${params.bbmap_ambiguous} \\
        minid=${params.bbmap_minid} \\
        pairlen=${params.bbmap_pairlen} \\
        path=index \\
        in=${reads[0]} \\
        in2=${reads[1]} \\
        out=${ref_id}_@_${pair_id}_bbmap_out.sam
    """
}

workflow {
    ref = channel.fromPath(params.test_bbmap_ref, checkIfExists: true)
        .map { ref_path ->
            def basename = ref_path.baseName
            def ref_id = basename.substring(0, basename.indexOf('.'))
            tuple(ref_id, ref_path)
        }

    reads = channel.fromFilePairs(params.test_bbmap_reads, checkIfExists: true)
   
    //ref.view { i -> "$i" }
    //reads.view { i -> "$i" }

    if (params.test_bbmap_separate) {
        index = bbmap_make_index(ref)
        //index.view { i -> "$i" }
        align_input = index.combine(reads)
        //align_input.view { i -> "$i" }
        bbmap_align = bbmap_align(align_input)
        //bbmap_align.view { i -> "$i" }
    } else {
        bbmap_input = ref.combine(reads)
        bbmap = bbmap(bbmap_input)
        //bbmap.view { i -> "$i" }
    }
}

