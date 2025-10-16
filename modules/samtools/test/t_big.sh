#!/bin/bash

n_threads=$(nproc)

for eid in 1620255 1620256 1620257 1620258
do
    sam=../../bwamem2/test/results.${eid}/out.sam

    wdir=work
    odir=results.${eid}
    log=t.log.${eid}

    mkdir -p ${odir} ${wdir}
    rm -rf ${odir}/* ${wdir}/* ${log}

    sam=$(cd $(dirname ${sam}) && pwd)/$(basename ${sam})
    wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

    apptainer exec --bind ${sam},${wdir},${odir} ../samtools.sif sh -c "\
              samtools view -@ ${n_threads} -b ${sam} > ${wdir}/unsorted.bam; \
              samtools sort -@ ${n_threads} -o ${odir}/out.bam ${wdir}/unsorted.bam; \
              samtools index -@ ${n_threads} ${odir}/out.bam; \
        " > ${log} 2>&1
done








