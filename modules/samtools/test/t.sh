#!/bin/bash -eu

#n_threads=$(nproc)
n_threads=16 # Since n_threads seems to affect the sorting result, it should be fixed here.

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

/usr/local/bin/apptainer exec --bind ${sam1},${odir} ../samtools.sif sh -c "\
    samtools sort -@ ${n_threads} -O bam -o ${odir}/out.sorted.bam ${sam1}" \
> ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${bam1},${wdir},${odir} ../samtools.sif sh -c "\
    samtools depth -m 10000000 -a ${bam1} -o ${wdir}/out.depth.txt \
    && head -n 100 ${wdir}/out.depth.txt \
    > ${odir}/out.depth.head100.txt" \
>> ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${bam1},${odir} ../samtools.sif sh -c "\
    samtools index ${bam1} -o ${odir}/out.index.bai" \
>> ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${bam1},${odir} ../samtools.sif sh -c "\
    samtools index ${bam1} -o ${wdir}/tmp.bai \
    && samtools bedcov ${bed1} ${bam1} -X ${wdir}/tmp.bai > ${odir}/out.bedcov.txt" \
>> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.txt ${refdir}/*.bam ${refdir}/*.bai)
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

# #{samtools_dir}/samtools sort -@ #{n_threads} -O bam -o #{out}.sort.bam #{out}.sam"
# #{samtools_dir}/samtools depth -m 10000000 -a #{out}.sort.bam
# #{samtools_dir}/samtools index #{read_cov}.sort.bam
# #{samtools_dir}/samtools bedcov #{orf_position}.bed #{read_cov}.sort.bam
# 10_contig_coverage.rb:samtools_dir









