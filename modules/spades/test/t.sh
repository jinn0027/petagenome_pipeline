#!/bin/bash

#n_threads=$(nproc)
n_threads=4
# Since #theads affects the results of spades, it's necessary to fix them.
mem=128

fq1=../../test/ecoli_1K_1.fq.gz
fq2=../../test/ecoli_1K_2.fq.gz

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py --only-error-correction --pe1-1 ${fq1} --pe1-2 ${fq2} --memory ${mem} --threads ${n_threads} -o ${odir} > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py --only-assembler --pe1-1 ${fq1} --pe1-2 ${fq2} --meta --memory ${mem} --threads ${n_threads} -o ${odir} >> ${log} 2>&1

for i in $(ls ${odir}/corrected/*.fastq.gz)
do
    gunzip $i
done

failed=""
for i in $(ls ${odir}/corrected/*.fastq)
do
    j=${refdir}/corrected/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
	failed="${failed} $(basename $i)"
        ret=1
    fi
done
for i in $(ls ${odir}/scaffolds.fasta)
do
    j=${refdir}/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
	failed="${failed} $(basename $i)"
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# #{spades_dir}/bin/spades.py --only-error-correction --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# #{spades_dir}/bin/spades.py --only-assembler --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --meta --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}
# 02_error_correction.rb 03_assembly.rb









