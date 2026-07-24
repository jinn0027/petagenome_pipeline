#!/bin/bash

#date=$(date +%Y%m%d)

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

PETAGENOME_PIPELINE_DIR=${MY_DIR}/..

memory=32
outdir=out
threads=16
cpus=16
random_seed=0

DIR_NF=$(pwd)/../nf
DIR_MODULES=$(pwd)/../modules
DIR_EXTERNAL=$(pwd)/../external

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

#args+=" --publish_output true"
#args+=" -profile slurm"

pushd ${DIR_EXTERNAL}

if [ ! -d bwa_db/GCF ] ; then
    if [ ! -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz ] ; then
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    fi

    nextflow run ${DIR_NF}/lv1/bwa_mem2.nf -entry BUILD_DB_ONLY ${args} \
	     --test_bwa_mem2_ref ${DIR_EXTERNAL}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

    mkdir -p bwa_db
    mv out/BUILD_DB_ONLY:BUILD_DB_SUB:bwa_mem2_makerefdb/00_GCF bwa_db/GRCh38
    rm -rf nfwork out
fi

popd
