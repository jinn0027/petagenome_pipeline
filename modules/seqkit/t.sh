#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} seqkit.sif seqkit -h

exit

# #{dir_tools_}/bin/seqkit replace --pattern "^(\\w+)(gene\\d+_gene_\\d+-)" --replacement "{kv}:${2}"
# #{seqkit_dir}/bin/seqkit grep -n --pattern-file <(cut -f 1 #{out_count}.txt) -o #{contig_proph}.fa
# 06_prophage_detection.02_after_processing.rb











