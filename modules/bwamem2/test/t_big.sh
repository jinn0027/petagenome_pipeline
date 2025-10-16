#!/bin/bash

n_threads=128

fq1=../../test/ecoli_1K_1.fq.gz
fq2=../../test/ecoli_1K_2.fq.gz
contig=../../megahit/test/ref/out/final.contigs.fa

for eid in 1620255 1620256 1620257 1620258
do
    fq1=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R1_001.fastq.gz
    fq2=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R2_001.fastq.gz
    contig=../../megahit/test/results.${eid}/out/final.contigs.fa

    wdir=work
    odir=results.${eid}

    fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
    fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
    wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

    log=t.log.${eid}

    mkdir -p ${wdir} ${odir}
    rm -rf ${wdir}/* ${odir}/* ${log}

    wref=${wdir}/$(basename ${contig})

    cp -f ${contig} ${wdir}

    date > ${log}
    echo "== Making index" >> ${log}
    apptainer exec --bind ${wdir} ../bwamem2.sif sh -c "\
              bwa-mem2 index ${wref}" >> ${log} 2>&1

    echo "== Alignment" >> ${log}
    date >> ${log}
    apptainer exec --bind ${fq1},${fq2},${wdir},${odir} ../bwamem2.sif sh -c "\
              bwa-mem2 mem -t ${n_threads} ${wref} ${fq1} ${fq2} > ${odir}/out.sam" >> ${log} 2>&1
    date >> ${log}

done
