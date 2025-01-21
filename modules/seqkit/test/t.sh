#!/bin/bash

fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../seqkit.sif seqkit -h > ${log} 2>&1
for i in $(ls $odir/*.html)
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

# #{dir_tools_}/bin/seqkit replace --pattern "^(\\w+)(gene\\d+_gene_\\d+-)" --replacement "{kv}:${2}"
# #{seqkit_dir}/bin/seqkit grep -n --pattern-file <(cut -f 1 #{out_count}.txt) -o #{contig_proph}.fa
# 06_prophage_detection.02_after_processing.rb











