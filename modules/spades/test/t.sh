#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} spades.sif /opt/spades/bin/spades.py --test --isolate

exit

# #{spades_dir}/bin/spades.py --only-error-correction --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# #{spades_dir}/bin/spades.py --only-assembler --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --meta --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# 02_error_correction.rb 03_assembly.rb









