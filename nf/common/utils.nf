#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def createNullParamsChannel() {
    return Channel.of(['@':'@'])
}

def getParam(obj, member) {
    return obj?.getAt(member) ?: params.getAt(member)
}

def clusterOptions(executor, gb, threads, label) {
    ret = ""
    envs = "PATH,LD_LIBRARY_PATH,PETAGENOME_PIPELINE_DIR"
    switch ("${executor}") {
        case "local" :
        break
    case "slurm" :
        ret += " --export=${envs}"
        if ("sc" in label) {
            ret +=" --exclusive"
        }
        break
    case "sge" :
        s_vmem = gb.toFloat() / threads.toFloat()
        ret += " -S /bin/bash -cwd -pe def_slot ${threads} -l s_vmem=${s_vmem}G -v ${envs}"
        break
    default :
        break
    }
    return ret
}

def processProfile(task) {
    return "### ${task.process} ${task.tag} on \${HOSTNAME}"
}

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

    def index_counter = new java.util.concurrent.atomic.AtomicInteger(0)
    def pairs = pairs_mixed.map { id, pair ->
        def index = index_counter.getAndIncrement()
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

    def index_counter = new java.util.concurrent.atomic.AtomicInteger(0)
    def seqs = seqs_mixed.map { id, seqs ->
        def index = index_counter.getAndIncrement()
        def new_id = "${String.format('%02d', index)}_${id}"
        return tuple(new_id, seqs)
    }

    return seqs
}
