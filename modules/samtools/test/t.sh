#!/bin/bash

# Since #threads may affect the sorting result, it should be fixed here.
#n_threads=$(nproc)
n_threads=16

sam1=../../test/minimiser-basic.sam
bam1=../../test/bedcov.bam
bed1=../../test/bedcov.bed

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

sam1=$(cd $(dirname ${sam1}) && pwd)/$(basename ${sam1})
bam1=$(cd $(dirname ${bam1}) && pwd)/$(basename ${bam1})
bed1=$(cd $(dirname ${bed1}) && pwd)/$(basename ${bed1})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
rm -rf ${odir}/* ${wdir}/*

apptainer exec --bind ${sam1},${wdir} ../samtools.sbx sh -c "\
    samtools sort -@ ${n_threads} -O bam -o ${wdir}/out.sorted.bam ${sam1}" \
> ${log} 2>&1

apptainer exec --bind ${sam1},${wdir},${odir} ../samtools.sbx sh -c "\
    samtools view -@ ${n_threads} -h ${wdir}/out.sorted.bam -o ${odir}/out.sorted.sam" \
> ${log} 2>&1

apptainer exec --bind ${bam1},${wdir},${odir} ../samtools.sbx sh -c "\
    samtools depth -m 10000000 -a ${bam1} -o ${wdir}/out.depth.txt \
    && head -n 100 ${wdir}/out.depth.txt \
    > ${odir}/out.depth.head100.txt" \
>> ${log} 2>&1

apptainer exec --bind ${bam1},${odir} ../samtools.sbx sh -c "\
    samtools index ${bam1} -o ${wdir}/tmp.bai \
    && samtools bedcov ${bed1} ${bam1} -X ${wdir}/tmp.bai > ${odir}/out.bedcov.txt" \
>> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.txt ${refdir}/*.sam)
do
    j=${odir}/$(basename $i)
    grep -v '@PG' $i > ${i}_
    grep -v '@PG' $j > ${j}_
    diff -q ${i}_ ${j}_ >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
	failed="${failed} $(basename $i)"
        ret=1
    fi
    rm ${i}_ ${j}_
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# #{samtools_dir}/samtools sort -@ #{n_threads} -O bam -o #{out}.sort.bam #{out}.sam"
# #{samtools_dir}/samtools depth -m 10000000 -a #{out}.sort.bam
# #{samtools_dir}/samtools index #{read_cov}.sort.bam
# #{samtools_dir}/samtools bedcov #{orf_position}.bed #{read_cov}.sort.bam
# 10_contig_coverage.rb:samtools_dir









