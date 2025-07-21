#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def createPairsChannel(path) {
    def pairs_list = path.split(';')

    def individual_channels = []
        pairs_list.each { pairs ->
        def ch = channel.fromFilePairs(pairs, checkIfExists: true)
        individual_channels << ch
    }

    def pairs_mixed = individual_channels.first()
    individual_channels.tail().each {
        ch -> pairs_mixed = pairs_mixed.mix(ch)
    }

    def index = 0
    def pairs = pairs_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index)}_${id}"
        index += 1
        return tuple(new_id, pair)
    }

    return pairs
}

def createSeqsChannel(path) {
    def seqs_list = path.split(';')

    def individual_channels = []
    seqs_list.each { seqs ->
        def ch = channel.fromPath(seqs, checkIfExists: true)
                        .collect()
                        .map{ it.sort() }
                        .map{ it ->
                            def key = it[0].simpleName.split('_')[0]
                            [key, it ]
                        }
        individual_channels << ch
    }

    def seqs_mixed = individual_channels.first()
    individual_channels.tail().each {
        ch -> seqs_mixed = seqs_mixed.mix(ch)
    }

    index = 0
    def seqs = seqs_mixed.map { id, seqs ->
        def new_id = "${String.format('%02d', index)}_${id}"
        index += 1
        return tuple(new_id, seqs)
    }

    return seqs
}