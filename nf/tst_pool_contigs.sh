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

inContigs=""
for i in ERR1620255 ERR1620256 ERR1620257 ERR1620258
do
    outDir=out_assembly_${i}/...
    if [ $i != ERR1620255 ] ; then
        inContigs+=" "
    fi
    inContigs+=""
done

nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
         -with-report report_${test}.html \
         -with-trace trace_${test}.txt \
         --test_pool_contigs_contigs "${inContigs}"
