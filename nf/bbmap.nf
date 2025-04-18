#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bbmap {
    container = "${params.petagenomeDir}/modules/bbmap/bbmap.sif"
    publishDir "${params.output}/bbmap/${ref_id}", mode: 'copy'
    input:
        tuple val(ref_id), path(ref)

    output:
        tuple val(ref_id), path("index")
    script:
    """
    bbmap.sh -Xmx${params.memory}g threads=${params.threads} ref=${ref} path=index
    """
}

workflow {
   test_bbmap_ref = channel.fromPath(params.test_bbmap_ref, checkIfExists: true)
        .map { ref_path ->
            def basename = ref_path.baseName
            def ref_id = basename.substring(0, basename.indexOf('.'))
            tuple(ref_id, ref_path)
        }

   test_bbmap_ref2 = test_bbmap_ref.map{a,b -> tuple(b,a)}

   //test_bbmap_reads = channel.fromFilePairs(params.test_bbmap_reads, checkIfExists: true)
   test_bbmap_ref.view { i -> "ref: ${i}" }
   test_bbmap_ref2.view { i -> "ref: ${i}" }
   //test_bbmap_reads.view { i -> "reads: ${i}" }
   //bbmap_out_index = bbmap(test_bbmap_ref)
}

