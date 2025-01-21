#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} prodigal.sif prodigal -h

exit

# #{dir_tools_}/Prodigal-GoogleImport/prodigal -a #{out_ex}.original.faa -d #{out_ex}.original.ffn -p meta -i #{contig_ex_fa_} -f gbk -o #{out_ex}.original.gbk
# 12_orf_prediction.from_05.rb 12_orf_prediction.from_06.rb







