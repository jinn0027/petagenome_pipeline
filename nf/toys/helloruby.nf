#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process helloruby {
    tag "${pair_id}"
    container = "${params.petagenomeDir}/modules/common/el9.sif"
    publishDir "${params.output}/helloruby/${pair_id}", mode: 'copy'
    // 当面は以下のように実行時にスクリプトパスを指定する。
    // 将来的にはel9.sifコンテナのビルド時にこのディレクトリごとバインドしてしまえばいいと思う。
    containerOptions = "--bind ${params.petagenomeDir}/scripts"
    input:
        tuple val(pair_id), path(reads, arity: '2')

    output:
        tuple val(pair_id), path("out/hello.txt")
    script:

        """
        mkdir -p out
        #ls /etc2 > out/hello.txt
        ruby ${params.petagenomeDir}/scripts/Ruby/toys/hello.rb > out/hello.txt
        """
}

workflow {
    reads = channel.fromFilePairs(params.test_helloruby_reads, checkIfExists: true)
    out = helloruby(reads)
    out.view { i -> "$i" }
}
