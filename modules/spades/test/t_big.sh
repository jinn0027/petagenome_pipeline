#!/bin/bash

n_threads=$(nproc)
mem=128

for eid in 1620255 1620256 1620257 1620258
do
    fq1=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R1_001.fastq.gz
    fq2=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R2_001.fastq.gz

    for k in erc asm all
    do
        odir=results.${eid}.${k}
        log=t.log.${eid}.${k}

        fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
        fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
        odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

        mkdir -p ${odir}
        rm -rf ${odir}/* ${log}

        date > ${log}
        if [ "$k" = "erc" ] ; then
            apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py \
                      --only-error-correction --pe1-1 ${fq1} --pe1-2 ${fq2} \
                      --memory ${mem} --threads ${n_threads} -o ${odir} > ${log} 2>&1
        elif [ "$k" = "asm" ] ; then
            apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py \
                      --only-assembler --pe1-1 ${fq1} --pe1-2 ${fq2} --meta \
                      --memory ${mem} --threads ${n_threads} -o ${odir} >> ${log} 2>&1
        else
            apptainer exec --bind ${fq1},${fq2},${odir} ../spades.sif spades.py \
                      --pe1-1 ${fq1} --pe1-2 ${fq2} --meta \
                      --memory ${mem} --threads ${n_threads} -o ${odir} >> ${log} 2>&1
        fi
        date >> ${log}
    done
done
