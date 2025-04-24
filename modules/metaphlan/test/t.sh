#!/bin/bash

# Since #threads may affect the result, it should be fixed here.
#n_threads=$(nproc)
n_threads=128

fq1=../../test/s_6_1.fastq.gz
#db=/opt/db

wdir=work
odir=results
refdir=ref
dbdir=../../../external/metaphlan4_db

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})

wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
dbdir=$(cd $(dirname ${dbdir}) && pwd)/$(basename ${dbdir})

mkdir -p ${odir} ${wdir}
rm -rf ${odir}/* ${wdir}/*

/usr/local/bin/apptainer exec --bind ${fq1},${odir},${dbdir} ../metaphlan.sif sh -c "\
    metaphlan \
        --bowtie2db ${dbdir} \
        --input_type fastq \
        --bowtie2out ${wdir}/out.all.sam \
        --nproc ${n_threads} \
        ${fq1} ${odir}/out.prof \
    && sort ${wdir}/out.all.sam | head -n 100 > ${odir}/out.sam \
    && sort ${wdir}/out.all.sam | tail -n 100 >> ${odir}/out.sam" > ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.sam)
do
    j=${odir}/$(basename $i)
    sort $i | grep -v '@PG' > ${wdir}/ii
    sort $j | grep -v '@PG' > ${wdir}/jj
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

# zcat #{query1_} #{query2_} | #{metaphlan_dir}/metaphlan2.py --bowtie2db ${DB} --input_type fastq --bowtie2out ${OUT_BWT2_} --nproc #{n_threads} > ${OUT_PROFILE_}
# 08_bacterial_taxonomic_profile.rb





