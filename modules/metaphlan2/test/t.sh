#!/bin/bash

n_threads=$(nproc)

fq1=../../test/s_6_1.fastq.gz
db=/opt/MetaPhlAn2/db_v20/mpa_v20_m200

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})

wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
rm -rf ${odir}/* ${wdir}/*

/usr/local/bin/apptainer exec --bind ${fq1},${odir} ../metaphlan2.sif sh -c "\
    zcat ${fq1} | \
    metaphlan2.py \
        --bowtie2db ${db} \
        --input_type fastq \
        --bowtie2out ${wdir}/out.all.sam \
        --nproc ${n_threads} \
        > ${odir}/out.prof \
    && head -n 100 ${wdir}/out.all.sam > ${odir}/out.sam \
    && tail -n 100 ${wdir}/out.all.sam >> ${odir}/out.sam" > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.prof ${refdir}/*.sam)
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

# zcat #{query1_} #{query2_} | #{metaphlan_dir}/metaphlan2.py --bowtie2db ${DB} --input_type fastq --bowtie2out ${OUT_BWT2_} --nproc #{n_threads} > ${OUT_PROFILE_}
# 08_bacterial_taxonomic_profile.rb











