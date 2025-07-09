#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512

test=pool_contigs

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"

for i in ERR1620255 ERR1620256 ERR1620257 ERR1620258
do
    ln -s out_assembly_${i}/assembly:filter_and_rename/out/contig.5000.fa contig.5000.${i}.fa
done

inContigs="contig.5000.*.fa"

outDir=out_pool_contigs

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outDir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    --publish_output true \
    "

nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
         -with-report report_${test}.html \
         -with-trace trace_${test}.txt \
         --test_pool_contigs_l_thre 5000 \
         --test_pool_contigs_contigs "${inContigs}"

rm -rf nfwork/*
