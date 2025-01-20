#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} samtools.sif samtools -h

exit

# #{samtools_dir}/samtools sort -@ #{n_threads} -O bam -o #{out}.sort.bam #{out}.sam
# #{samtools_dir}/samtools depth -m 10000000 -a #{out}.sort.bam
# #{samtools_dir}/samtools index #{read_cov}.sort.bam
# #{samtools_dir}/samtools bedcov #{orf_position}.bed #{read_cov}.sort.bam
# 10_contig_coverage.rb:samtools_dir









