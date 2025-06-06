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
rm -rf ${odir}/*

apptainer exec --bind ${fq1},${fq2},${odir} ../fastqc.sbx fastqc -o ${odir} ${fq1} ${fq2} > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.html)
do
    j=${odir}/$(basename $i)
    sed 's#header_filename">[^<]*<##g' $i > _ref
    sed 's#header_filename">[^<]*<##g' $j > _out
    diff -q _ref _out >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        failed="${failed} $(basename $i)"
        ret=1
    fi
    rm -f _ref _out
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# #{fastqc_dir}/bin/fastqc -o #{dir_fastqc_} #{out_cutadapt}_1.fastq
# 01_trim_qc.HiSeq.rb 02_error_correction.rb

