#!/bin/bash

if [[ ! -v PETAGENOME_PIPELINE_DIR ]] ; then
    echo "PETAGENOME_PIPELINE_DIR not defined"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
if [ ! -d ${PETAGENOME_PIPELINE_DIR} ] ; then
    echo "${PETAGENOME_PIPELINE_DIR} does not exist"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
echo "PETAGENOME_PIPELINE_DIR : ${PETAGENOME_PIPELINE_DIR}"

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

date=$(date +"%Y%m%d%H%M%S")

threads=1 #$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=10

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
outDir=out
inPairs="out/fastp/*/out_{1,2}.fastq"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outDir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    --publish_output true \
    "

args+="\
    --fastp_fastp_memory 10 \
    --fastp_fastp_threads 30 \
    --spades_spades_error_correction_memory 80 \
    --spades_spades_error_correction_threads 30 \
    --spades_spades_assembler_memory 100 \
    --spades_spades_assembler_threads 30 \
    --mmseqs2_mmseqs2_cluster_memory 100 \
    --mmseqs2_mmseqs2_cluster_threads 100 \
    "

args_dbg="\
    -with-report report_error_correction.${date}.html \
    -with-trace trace_error_correction.${date}.txt \
    -with-timeline timeline_error_correction.${date}.html \
    -with-dag dag_error_correction.${date}.png \
    "

nextflow clean -f
nextflow run ${nfDir}/lv2/error_correction.nf \
    ${args} \
    ${args_dbg} \
    --test_error_correction_reads "${inPairs}"

#rm -rf nfwork/*
