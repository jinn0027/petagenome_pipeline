#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

date=$(date +"%Y%m%d%H%M%S")

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512

lthre=5000

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
inPairs="/scratch/local/data/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"

outDir=out

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
         -with-report report_bacteriome_pipeline.${date}.html \
         -with-trace trace_bacteriome_pipeline.${date}.txt \
         -with-timeline timeline_bacteriome_pipeline.${date}.html \
         -with-dag dag_bacteriome_pipeline.${date}.png \
         --test_bacteriome_pipeline_lthre "${lthre}" \
         --test_bacteriome_pipeline_reads "${inPairs}"

#rm -rf nfwork/*

