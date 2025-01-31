#!/bin/bash

# Since #threads may affect the sorting result, it should be fixed here.
#n_threads=$(nproc)
n_threads=128
random_seed=1

fa1=../../test/1.faa

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
rm -rf ${odir}/* ${wdir}/*

pfam_hmm=/opt/Pfam-A.hmm

/usr/local/bin/apptainer exec --bind ${fa1},${odir} ../hmmer.sbx sh -c "\
    hmmstat ${pfam_hmm} > ${odir}/pfam_stat.txt" > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${fa1},${odir} ../hmmer.sbx sh -c "\
    hmmscan --seed ${random_seed} --cpu ${n_threads} --domtblout ${odir}/sample_hmmscan.txt ${pfam_hmm} ${fa1}" > ${log} 2>&1

failed=""

for i in $(ls ${refdir}/*.txt)
do
    j=${odir}/$(basename $i)
    grep -v '#' $i > ${wdir}/ii
    grep -v '#' $j > ${wdir}/jj
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

# NOP : results seem to be used but command is not found; seems to be used in virsorter etc











