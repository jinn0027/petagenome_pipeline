#!/bin/bash

vcf1=../../test/idx.vcf

odir=results
refdir=ref

log=t.log

ret=0

vcf1=$(cd $(dirname ${vcf1}) && pwd)/$(basename ${vcf1})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

# get statics
/usr/local/bin/apptainer exec --bind ${vcf1},${odir} ../bcftools.sif sh -c "\
    bcftools stats ${vcf1} > ${odir}/stats.txt" > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.txt)
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

# NOP












