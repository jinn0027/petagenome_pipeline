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
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../bedtools.sif bedtools -h  2>&1
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

# #{bedtools_dir}/bin/bedtools merge -i #{blast_out}.bed -d 150 -s -c 4,5,6,7,8,9 -o distinct,distinct,distinct,distinct,min,max
# #{bedtools_dir}/bin/bedtools merge -i #{blast_out}.merged.tmp2.bed -d 0 -c 4,5,6 -o distinct,distinct,distinct
# #{bedtools_dir}/bin/bedtools getfasta -fi #{contig_proph}.fa -bed #{blast_out}.merged.bed -fo #{merged_prophage}.original.fa
# #{bedtools_dir}/bin/bedtools intersect -a #{orf_position}.bed -b #{read_cov}.region.bed -wo > #{out_orf}.intersect.bed
# 06_prophage_detection.02_after_processing.rb 10_orf_coverage.rb











