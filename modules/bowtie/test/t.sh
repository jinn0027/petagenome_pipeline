#!/bin/bash

#n_threads=$(nproc)
n_threads=4
# Since #theads may affect the results of bowtie, it's necessary to fix them.
random_seed=1

fa1=../../test/8seq.fa
fa2=../../test/1seq.fa

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
fa2=$(cd $(dirname ${fa2}) && pwd)/$(basename ${fa2})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

idx=${wdir}/$(basename ${fa1})

mkdir -p ${wdir} ${odir}
rm -rf ${wdir}/* ${odir}/*

# create index
apptainer exec --bind ${fa1},${wdir} ../bowtie.sbx sh -c "\
    bowtie-build ${fa1} ${idx} --threads ${n_threads} --seed ${random_seed}" > ${log} 2>&1

# align
apptainer exec --bind ${fa2},${wdir},${odir} ../bowtie.sbx sh -c "\
    bowtie -S -f -x ${idx} ${fa2} -p ${n_threads} --seed ${random_seed} > ${odir}/out.sam" >> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.sam)
do
    j=${odir}/$(basename $i)
    head -n 10 $i | grep -v '@PG' > ${wdir}/ii
    head -n 10 $j | grep -v '@PG' > ${wdir}/jj
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

# NOP : results seem to be used but command is not found








