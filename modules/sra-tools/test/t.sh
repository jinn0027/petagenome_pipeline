#!/bin/bash

odir=results
refdir=ref

log=t.log

ret=0

odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

/usr/local/bin/apptainer exec --bind ${odir} ../sra-tools.sbx sh -c "\
    fastq-dump --stdout -X 2 SRR390728 > ${odir}/SRR390728.txt" > ${log} 2>&1

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










