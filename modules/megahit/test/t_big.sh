#!/bin/bash

n_threads=$(nproc)
#n_threads=4
# Since #theads may affect the results of megahit, it's necessary to fix them.
mem=128

#fq1=../../test/ecoli_1K_1.fq.gz
#fq2=../../test/ecoli_1K_2.fq.gz
fq1=./ERR1620255_XXXXXXXX_XXXXXXXX_L001_R1_001.fastq.gz
fq2=./ERR1620255_XXXXXXXX_XXXXXXXX_L001_R2_001.fastq.gz

odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

mem_in_byte=${mem}000000000

apptainer exec --bind ${fq1},${fq2},${odir} ../megahit.sif \
    megahit -1 ${fq1} -2 ${fq2} -o ${odir}/out -m ${mem_in_byte} -t ${n_threads} > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/out/*.fa)
do
    j=${odir}/out/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        failed="${failed} $(basename $i)"
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# NOP











