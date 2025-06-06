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

apptainer exec --bind ${fa1},${fa2},${wdir} ../mmseqs2.sbx sh -c "\
    mmseqs createdb ${fa1} ${wdir}/$(basename ${fa1})" > ${log} 2>&1

apptainer exec --bind ${fa1},${fa2},${wdir} ../mmseqs2.sbx sh -c "\
    mmseqs createdb ${fa2} ${wdir}/$(basename ${fa2})" > ${log} 2>&1

apptainer exec --bind ${fa1},${fa2},${wdir},${odir} ../mmseqs2.sbx sh -c "\
    mmseqs search --search-type 3 --threads ${n_threads} ${wdir}/$(basename ${fa2}) ${wdir}/$(basename ${fa1}) ${odir}/out ${wdir}" > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*)
do
    ii=$(basename ${i})
    diff -q ${odir}/${ii} ${refdir}/${ii} >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        failed="${failed} ${ii}"
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










