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
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../samtools.sif samtools -h > ${log} 2>&1
for i in $(ls $odir/*.html)
do
    j=$refdir/$(basename $i)
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

# #{samtools_dir}/samtools sort -@ #{n_threads} -O bam -o #{out}.sort.bam #{out}.sam
# #{samtools_dir}/samtools depth -m 10000000 -a #{out}.sort.bam
# #{samtools_dir}/samtools index #{read_cov}.sort.bam
# #{samtools_dir}/samtools bedcov #{orf_position}.bed #{read_cov}.sort.bam
# 10_contig_coverage.rb:samtools_dir









