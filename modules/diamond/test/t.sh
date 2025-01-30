#!/bin/bash

n_threads=$(nproc)

fa1=../../test/1.faa
fa2=../../test/2.faa

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
fa2=$(cd $(dirname ${fa2}) && pwd)/$(basename ${fa2})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

db=${wdir}/db

mkdir -p ${wdir} ${odir}
rm -rf ${wdir}/* ${odir}/*

/usr/local/bin/apptainer exec --bind ${fa1},${fa2},${wdir} ../diamond.sbx sh -c "\
  diamond makedb --threads ${n_threads} --in ${fa1} -d ${db}" > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${fa1},${fa2},${wdir} ../diamond.sbx sh -c "\
  diamond blastp -d ${db} -q ${fa2} -o ${odir}/out.tsv" >> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.tsv)
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










