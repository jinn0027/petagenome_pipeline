#!/bin/bash

fa1=../../test/ecoli_1K_1.fa.gz

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}

wfa1=${wdir}/$(basename ${fa1} | sed 's#.gz$##')

gunzip ${fa1} -c > ${wfa1}

/usr/local/bin/apptainer exec --bind ${wfa1},${odir} ../prodigal.sif prodigal -a ${odir}/out.faa -d ${odir}/out.ffn -p meta -i ${wfa1} -f gbk -o ${odir}/out.gbk > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.faa ${refdir}/*.ffn ${refdir}/*.gbk)
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

# #{dir_tools_}/Prodigal-GoogleImport/prodigal -a #{out_ex}.original.faa -d #{out_ex}.original.ffn -p meta -i #{contig_ex_fa_} -f gbk -o #{out_ex}.original.gbk
# 12_orf_prediction.from_05.rb 12_orf_prediction.from_06.rb







