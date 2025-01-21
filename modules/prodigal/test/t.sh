#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../prodigal.sif prodigal -h > ${log} 2>&1

failed=""
for i in $(ls ${odir}/*.fa)
do
    j=${refdir}/$(basename $i)
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







