#!/bin/bash

# Since #threads may affect the result, it should be fixed here.
#n_threads=$(nproc)
n_threads=128

fq1=../../test/ecoli_1K_1.fq.gz
fq2=../../test/ecoli_1K_2.fq.gz
contig=../../megahit/test/ref/out/final.contigs.fa
config_templ=./config.templ

odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
contig=$(cd $(dirname ${contig}) && pwd)/$(basename ${contig})
config_templ=$(cd $(dirname ${config_templ}) && pwd)/$(basename ${config_templ})

odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

apptainer exec --bind ${fq1},${fq2},${contig},${config_templ},${odir} ../soapdenovo2.sif sh -c "\
  cat ${config_templ} > config ; \
  echo \"q1=${fq1}\" >> config ; \
  echo \"q2=${fq2}\" >> config ; \
  cat config ; \
  SOAPdenovo-fusion -D -s config -p 40 -K 41 -g k41 -c ${contig}; \
  SOAPdenovo-127mer map -s config -p 40 -g k41; \
  SOAPdenovo-127mer scaff -p 40 -g k41; \
  ls; \
"

exit

  cp ${config_templ} ./config \

##-t ${n_threads} -o ${odir} ${fq1} ${fq2} > ${log} 2>&1

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










