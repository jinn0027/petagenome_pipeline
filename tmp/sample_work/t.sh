#!/bin/bash

export TMPDIR=/dev/shm/tmp

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=128 # GB
outdir=/dev/shm/${USER}/petagenome_pipeline/out

petagenomeDir=/scratch/local/petagenome_pipeline

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

