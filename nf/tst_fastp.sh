#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512

test=fastp

outdir=out_${test}

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"
extDir="${PETAGENOME_PIPELINE_DIR}/external"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

args+=" --publish_output true"

dataDir=/scratch/local/data/metagenome

for i in ERR1620255 ERR1620256 ERR1620257 ERR1620258
do
    inPairs="${dataDir}/${i}_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
    nextflow run ${nfDir}/lv1/fastp.nf ${args} \
             -with-report report_${test}_${i}.html \
             -with-trace trace_${test}_${i}.txt \
             --test_fastp_reads "${inPairs}"
done
