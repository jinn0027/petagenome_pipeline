#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} bowtie2.sif bowtie2 -h

exit

# NOP : results seem to be used but command is not found








