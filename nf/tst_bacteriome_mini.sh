#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

date=$(date +"%Y%m%d%H%M%S")

threads=10 #$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=10

#lthre=5000
lthre=0

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
#inPairs="/scratch/local/data/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
inPairs="${PETAGENOME_PIPELINE_DIR}/modules/test/ecoli_1K_{1,2}.fq.gz:${PETAGENOME_PIPELINE_DIR}/modules/test/s_6_{1,2}.fastq.gz"

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

#args+="\
#    --spades_spades_error_correction_memory 80 \
#    --spades_spades_error_correction_threads 30 \
#    --spades_spades_assembler_memory 100 \
#    --spades_spades_assembler_threads 30 \
#    --mmseqs2_mmseqs2_cluster_memory 100 \
#    --mmseqs2_mmseqs2_cluster_threads 100 \
#    "

nextflow clean -f
nextflow run ${nfDir}/lv3/bacteriome_pipeline.nf ${args} \
         -with-trace trace_bacteriome_pipeline.${date}.txt \
         --test_bacteriome_pipeline_lthre "${lthre}" \
         --test_bacteriome_pipeline_reads "${inPairs}"

#         -with-report report_bacteriome_pipeline.${date}.html \

#rm -rf nfwork/*

