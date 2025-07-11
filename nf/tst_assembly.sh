#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512

test=assembly

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"

for i in ERR1620255 ERR1620256 ERR1620257 ERR1620258
do
    inDir=out_error_correction_${i}/error_correction:spades_error_correction/out/corrected/paired/
    inPairs="${inDir}/out_{1,2}00.0_0.cor.fastq"
    outDir=out_assembly_${i}

    args="\
        --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
        --output ${outDir} \
        --memory ${memory} \
        --threads ${threads} \
        --cpus ${cpus} \
        --random_seed ${random_seed} \
        --publish_output true \
        "

    nextflow run ${nfDir}/lv2/assembly.nf ${args} \
             -with-report report_${test}_${i}.html \
             -with-trace trace_${test}_${i}.txt \
             --test_assembly_l_thre 5000 \
             --test_assembly_reads "${inPairs}"

    rm -rf nfwork/*
done
