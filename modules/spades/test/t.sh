#!/bin/bash

n_threads=$(nproc)
mem=128

fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py --only-error-correction --pe1-1 ${fq1} --pe1-2 ${fq2} --memory ${mem} --threads ${n_threads} -o ${odir} > ${log} 2>&1
for i in $(ls ${odir}/corrected "*.fastq.gz")
do
    j=${refdir}/corrected/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED" | tee -a ${log}
fi

exit ${ret}

# #{spades_dir}/bin/spades.py --only-error-correction --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# #{spades_dir}/bin/spades.py --only-assembler --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --meta --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# 02_error_correction.rb 03_assembly.rb









