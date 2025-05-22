#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

#threads=$(nproc)
threads=16
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
#memory=128 # GB
memory=32 # GB
#outdir=/dev/shm/${USER}/petagenome_pipeline/out
outdir=out

if [ -v PETAGENOME_PIPELINE_DIR ] && [ -d ${PETAGENOME_PIPELINE_DIR} ] ; then
    petagenomeDir=${PETAGENOME_PIPELINE_DIR}
elif [ -v HOME ] && [ -d ${HOME}/petagenome_pipeline ] ; then
    petagenomeDir=${HOME}/petagenome_pipeline
elif [ -d /scratch/local/petagenome_pipeline ] ; then
    petagenomeDir=/scratch/local/petagenome_pipeline
fi

if [ "${petagenomeDir}" = "" ] ; then
    echo "petagenome_pipeline not found"
    exit 1
fi

args="\
    --petagenomeDir=${petagenomeDir} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

# ローカルのmain.nfを実行
nextflow run main.nf ${args} --test_main_reads 'ecoli_1K_{1,2}.fq.gz'

# 以下は/scratch/local/petagenome_pipeline/nf以下のNextflowスクリプトを実行する場合

#nextflow run ${petagenomeDir}/nf/main.nf ${args} --test_main_reads 'ecoli_1K_{1,2}.fq.gz'
#nextflow run ${petagenomeDir}/nf/cutadapt.nf ${args} --test_cutadapt_reads 'ecoli_1K_{1,2}.fq.gz'

