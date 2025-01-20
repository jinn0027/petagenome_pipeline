#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} virsorter.sif virsorter -h

exit

# perl #{virsorter_dir}/wrapper_phage_contigs_sorter_iPlant.pl #{opt_virome} --db #{vs_db_id} -data-dir #{dir_data_vs_} -ncpu #{n_threads} -wdir #{dir_out_} -f #{query_fa_}
# 06_prophage_detection.01_each_sample.rb












