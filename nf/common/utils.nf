#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def createReadsChannel(reads_param) {
    def reads_list = reads_param.split(';')

    def individual_channels = []
        reads_list.each { reads ->
        def ch = channel.fromFilePairs(reads, checkIfExists: true)
        individual_channels << ch
    }

    def reads_mixed = individual_channels.first()
    individual_channels.tail().each {
        ch -> reads_mixed = reads_mixed.mix(ch)
    }

    def index = 0
    def reads = reads_mixed.map { id, pair ->
        def new_id = "${String.format('%02d', index)}_${id}"
        index += 1
        return tuple(new_id, pair)
    }

    return reads
}

def createPathChannel(path_param) {
    def path_list = path_param.split(';')

    def individual_channels = []
    path_list.each { path ->
        def ch = channel.fromPath(path, checkIfExists: true)
                    .collect()
                    .map{ it.sort() }
                    .map{ it ->
                          def key = it[0].simpleName.split('_')[0]
                     [key, it ] }
        individual_channels << ch
    }

    def path_mixed = individual_channels.first()
    individual_channels.tail().each {
        ch -> path_mixed = path_mixed.mix(ch)
    }

    index = 0
    def path = path_mixed.map { id, path ->
        def new_id = "${String.format('%02d', index)}_${id}"
        index += 1
        return tuple(new_id, path)
    }

    return path
}