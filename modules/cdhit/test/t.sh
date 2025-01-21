#!/bin/bash

n_threads=$(nproc)
mem=128

fa1=../../test/s_6_1.fasta.gz
fa2=../../test/s_6_2.fasta.gz
odir=results

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
fa2=$(cd $(dirname ${fa2}) && pwd)/$(basename ${fa2})

odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
refdir=ref

log=t.log

ofa1=${odir}/$(basename ${fa1} | sed 's#.gz$##')
ofa2=${odir}/$(basename ${fa2} | sed 's#.gz$##')

ret=0

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fa1},${fa2},${odir} ../cdhit.sif cd-hit-est -c 0.95 -G 1 -mask NX -d 150 -n 10 -T ${n_threads} -M ${mem}000 -i ${fa1} -o ${ofa1} > ${log} 2>&1
for i in $(ls $odir/*.fasta)
do
    j=${refdir}/$(basename $i)
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

# #{cdhit_dir}/cd-hit-est -c 0.95 -G 1 -mask NX -d 150 -n 10 -T #{n_threads} -M #{mem}000 -i #{merged_fa_} -o #{rep_fa_}
# 04_pool_contigs.rb











