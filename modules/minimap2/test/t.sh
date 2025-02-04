#!/bin/bash

# Since #threads may affect the result, it should be fixed here.
#n_threads=$(nproc)
n_threads=128

fa1=../../test/8seq.fa
fa2=../../test/1seq.fa

wdir=work
odir=results

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
fa2=$(cd $(dirname ${fa2}) && pwd)/$(basename ${fa2})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
refdir=ref

log=t.log

ret=0

mkdir -p ${wdir} ${odir}
rm -rf ${wdir}/* ${odir}/*

/usr/local/bin/apptainer exec --bind ${fa1},${fa2},${odir} ../minimap2.sbx sh -c "\
    minimap2 -t ${n_threads} -a ${fa1} ${fa2} > ${odir}/out.sam" > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.sam)
do
    j=${odir}/$(basename $i)
    sort $i | grep -v '@PG' > ${wdir}/ii
    sort $j | grep -v '@PG' > ${wdir}/jj
    diff -q ${wdir}/ii ${wdir}/jj >> ${log} 2>&1 && :
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










