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

mkdir -p uniprot_refs
pushd uniprot_refs
if [ ! -f uniprot_sprot.fasta.gz  ] ; then
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz 
fi

if [ ! -f uniprot_to_taxid.tsv  ] || [ ! -f uniprot_to_ko.tsv ] ; then
    if [ ! -f idmapping_selected.tab.gz ] ; then
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    fi
    zcat idmapping_selected.tab.gz | awk -F'\t' ' { print $1 "\t" $12 > "uniprot_to_ko.tsv"; print $1 "\t" $3 > "uniprot_to_taxid.tsv"; } '
    sed -i 's#; #;#g' uniprot_to_ko.tsv
    sed -i 's#; #;#g' uniprot_to_taxid.tsv
fi
popd
popd
