#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} cdhit.sif cdhit -h

exit

# #{cdhit_dir}/cd-hit-est -c 0.95 -G 1 -mask NX -d 150 -n 10 -T #{n_threads} -M #{mem}000 -i #{merged_fa_} -o #{rep_fa_}
# 04_pool_contigs.rb











