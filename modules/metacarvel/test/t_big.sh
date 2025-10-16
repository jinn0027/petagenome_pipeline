#!/bin/bash

n_threads=128


for eid in 1620255 1620256 1620257 1620258
do
    bam=../../samtools/test/results.${eid}/out.bam
    contig=../../megahit/test/results.${eid}/out/final.contigs.fa

    odir=results.${eid}
    log=t.log.${eid}

    mkdir -p ${odir}
    rm -rf ${odir}/* ${log}

    bam=$(cd $(dirname ${bam}) && pwd)/$(basename ${bam})
    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
    contig=$(cd $(dirname ${contig}) && pwd)/$(basename ${contig})

    date > ${log}
    apptainer exec --bind ${contig},${bam},${odir} \
        ../metacarvel.sif sh -c "\
            python /opt/MetaCarvel/run.py \
                -a ${contig} \
                -m ${bam} \
                -d ${odir} ; \
        " >> ${log} 2>&1
    date >> ${log}
done


