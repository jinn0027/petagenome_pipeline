#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

#threads=16
threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=512
#memory=128
#outdir=/dev/shm/${USER}/petagenome_pipeline/out
outdir=out

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"
extDir="${PETAGENOME_PIPELINE_DIR}/external"

sample="ERR1620255"
#sample="ERR1620256"
#sample="ERR1620257"
#sample="ERR1620258"

testPair="/scratch/local/data/metagenome/${sample}_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
l_thre=5000
fastpOutPair="out/fastp/ERR1620255_XXXXXXXX_XXXXXXXX_L001_R/out_{1,2}.fastq*"
errorCorrectionOutPair="out/error_correction:spades_error_correction/out/corrected/paired/out_{1,2}00.0_0.cor.fastq"
assemblyOutContigs="out/assembly_filter_and_rename/out/contig.${l_thre}.fa"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

args+=" --publish_output true"

#test=${1:-"fastp"}
#test=${1:-"error_correction"}
#test=${1:-"assembly"}
test=${1:-"pool_contigs"}
##test=${1:-"circular_contigs"}

test=${test%.*}

case ${test} in
    "main")
        nextflow run ${nfDir}/toys/main.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_main_reads "${testPair}"
        ;;
    "fastp")
        nextflow run ${nfDir}/lv1/fastp.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_fastp_reads "${testPair}"
        ;;
    "error_correction")
        nextflow run ${nfDir}/lv2/error_correction.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_error_correction_reads "${fastpOutPair}"
        ;;
    "assembly")
        nextflow run ${nfDir}/lv2/assembly.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_assembly_reads "${errorCorrectionOutPair}"
        ;;
    "pool_contigs")
        nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_pool_contigs_contigs "${assemblyOutContigs}"
        ;;
    "*")
esac

#rm -rf /dev/shm/${USER}
