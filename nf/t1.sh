#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=16
#threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=256
#memory=128
#outdir=/dev/shm/${USER}/petagenome_pipeline/out
outdir=out

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"
extDir="${PETAGENOME_PIPELINE_DIR}/external"

testPair1="/scratch/local/data/metagenome/ERR1620255_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
testPair2="/scratch/local/data/metagenome/ERR1620256_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
testPair3="/scratch/local/data/metagenome/ERR1620257_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
testPair4="/scratch/local/data/metagenome/ERR1620258_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

args+=" --publish_output true"

xargs+="\
    -with-trace \
    -with-report \
    "

#test=${1:-"fastp"}
#test=${1:-"error_correction"}
#test=${1:-"assembly"}
test=${1:-"pool_contigs"}
##test=${1:-"circular_contigs"}

test=${test%.*}

case ${test} in
    "main")
        nextflow run ${nfDir}/toys/main.nf ${args} \
                 --test_main_reads "${testPair1}"
        ;;
    "fastp")
        nextflow run ${nfDir}/lv1/fastp.nf ${args} \
                 --test_fastp_reads "${testPair1}"
        ;;
    "error_correction")
        nextflow run ${nfDir}/lv2/error_correction.nf ${args} \
                 --test_error_correction_reads "${testPair1}"
        ;;
    "assembly")
        nextflow run ${nfDir}/lv2/assembly.nf ${args} \
                 --test_assembly_reads "${testPair1}"
        ;;
    "pool_contigs")
        nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
                 --test_pool_contigs_contigs "${testPair1}"
        ;;
    "*")
esac

rm -rf /dev/shm/${USER}

#                 -with-report report.html \
#                 -with-trace \
