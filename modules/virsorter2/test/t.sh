#!/bin/bash

# Since #threads may affect the result, it should be fixed here.
#n_threads=$(nproc)
n_threads=128

fa1=../../test/1seq.fa

odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

/usr/local/bin/apptainer exec --bind ${fa1},${odir} --writable ../virsorter2.sbx virsorter run -w ${odir} -i ${fa1} -j ${n_threads}  all > ${log} 2>&1

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







