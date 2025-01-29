#!/bin/bash

fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

odir=results
refdir=ref

log=t.log

ofq1=${odir}/$(basename ${fq1} | sed 's#.gz$##')
ofq2=${odir}/$(basename ${fq2} | sed 's#.gz$##')

ret=0

ADPT_FWD="AATGATACGGCGACCACCGAGAUCTACAC"
ADPT_REV="CAAGCAGAAGACGGCATACGAGAT"

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../cutadapt.sif cutadapt --minimum-length 50 -a ${ADPT_FWD} -g ${ADPT_REV} -o ${ofq1} -p ${ofq2} ${fq1} ${fq2}> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.fastq)
do
    j=${odir}/$(basename $i)
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

# #{cutadapt_dir}/bin/cutadapt --minimum-length 50 -a ${ADPT_FWD} -o #{out_cutadapt}_1.fastq ${QUERY1_}
# 01_trim_qc.HiSeq.rb




