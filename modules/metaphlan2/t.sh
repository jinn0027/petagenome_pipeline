#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} metaphlan2.sif metaphlan2.py --help

exit

# #{metaphlan_dir}/metaphlan2.py --bowtie2db ${DB} --input_type fastq --bowtie2out ${OUT_BWT2_} --nproc #{n_threads}
# 08_bacterial_taxonomic_profile.rb











