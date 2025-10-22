#!/bin/bash

petagenomeDir=""
if [[ -v PETAGENOME_PIPELINE_DIR ]] ; then
    petagenomeDir=${PETAGENOME_PIPELINE_DIR}
else
    echo "PETAGENOME_PIPELINE_DIR not defined"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
if [ ! -d ${petagenomeDir} ] ; then
    echo "${petagenomeDir} does not exist"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
echo "PETAGENOME_PIPELINE_DIR : ${petagenomeDir}"

pzrepoDir=""
if [[ -v PZREPO_DIR ]] ; then
    pzrepoDir=${PZREPO_DIR}
else
    echo "PZREPO_DIR not defined"
    echo "If necessary, please source <pzrepo_dir>/etc/host_setup.sh"
fi
if [ ! -d ${pzrepoDir} ] ; then
    echo "${pzrepoDir} does not exist"
    echo "Please source <pzrepo_dir>/etc/host_setup.sh"
    exit 1
fi
echo "PZREPO_DIR : ${pzrepoDir}"

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

mkdir -p ${TMPDIR}

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

date=$(date +"%Y%m%d%H%M%S")

threads=10 #$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=20

lthre=5000

nfDir="${petagenomeDir}/nf"
inPairs="/scratch/local/data/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"
#inPairs="~/data/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"

outDir=out

args="\
    --petagenomeDir=${petagenomeDir} \
    --output ${outDir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    --publish_output true \
    "

if [ "${pzrepoDir}" != "" ] ; then
    args+=" --pzrepoDir=${pzrepoDir}"
fi

args+="\
    -profile local
    "

args+="\
    --fastp_fastp_memory 30 \
    --fastp_fastp_threads 30 \
    --spades_spades_error_correction_memory 80 \
    --spades_spades_error_correction_threads 30 \
    --spades_spades_assembler_memory 100 \
    --spades_spades_assembler_threads 30 \
    --mmseqs2_mmseqs2_cluster_memory 100 \
    --mmseqs2_mmseqs2_cluster_threads 100 \
    "

args_dbg="\
    -with-report report_bacteriome_pipeline.${date}.html \
    -with-trace trace_bacteriome_pipeline.${date}.txt \
    -with-timeline timeline_bacteriome_pipeline.${date}.html \
    -with-dag dag_bacteriome_pipeline.${date}.png \
    "

nextflow clean -f
nextflow run ${nfDir}/lv3/bacteriome_pipeline.nf \
    ${args} \
    ${args_dbg} \
    --test_bacteriome_pipeline_lthre "${lthre}" \
    --test_bacteriome_pipeline_reads "${inPairs}"

#rm -rf nfwork/*

