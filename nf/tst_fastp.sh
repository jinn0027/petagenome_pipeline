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

threads=1 # $(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=10

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
outDir=out
inPairs="/scratch/local/data/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
#inPairs="${PETAGENOME_PIPELINE_DIR}/modules/test/ecoli_1K_{1,2}.fq.gz;${PETAGENOME_PIPELINE_DIR}/modules/test/s_6_{1,2}.fastq.gz"

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
    -with-report report_fastp.${date}.html \
    -with-trace trace_fastp.${date}.txt \
    -with-timeline timeline_fastp.${date}.html \
    -with-dag dag_fastp.${date}.png \
    "

nextflow clean -f
nextflow run ${nfDir}/lv1/fastp.nf \
    ${args} \
    ${args_dbg} \
    --test_fastp_reads "${inPairs}"

#rm -rf nfwork/*

