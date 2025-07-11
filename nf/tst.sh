#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512

test=bacteriome_pipeline

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"

#dataDir=
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"

inPairs="/scratch/local/data/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
#inPairs="${PETAGENOME_PIPELINE_DIR}/modules/test/ecoli_1K_{1,2}.fq.gz"

outDir=out_bacteriome_pipeline

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outDir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    --publish_output true \
    "

nextflow clean -f
nextflow run ${nfDir}/lv3/bacteriome_pipeline.nf ${args} \
         -with-report report_${test}.html \
         -with-trace trace_${test}.txt \
         --test_bacteriome_pipeline_reads "${inPairs}"

#rm -rf nfwork/*

