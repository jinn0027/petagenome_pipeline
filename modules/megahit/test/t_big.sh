#!/bin/bash

n_threads=$(nproc)
mem=128
mem_in_byte=${mem}000000000

for eid in 1620255 1620256 1620257 1620258
do
    fq1=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R1_001.fastq.gz
    fq2=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R2_001.fastq.gz

    odir=results.${eid}
    log=t.log.${eid}

    mkdir -p ${odir}
    rm -rf ${odir}/* ${log}

    fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
    fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

    date > ${log}
    apptainer exec --bind ${fq1},${fq2},${odir} ../megahit.sif \
              megahit -1 ${fq1} -2 ${fq2} -o ${odir}/out -m ${mem_in_byte} -t ${n_threads} >> ${log} 2>&1
    date >> ${log}
done

