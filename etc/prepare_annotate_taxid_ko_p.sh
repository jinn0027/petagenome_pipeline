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

# -----------------------------------------------------------------
# 1. KO ID -> KO Definition (KEGG API から取得)
# -----------------------------------------------------------------
if [ ! -f ko_to_name.tsv ]; then
    echo "Downloading KEGG KO list..."
    # KEGG REST API から ko と Definition の対応表を取得
    # 出力形式: ko:K00001\talcohol dehydrogenase [EC:1.1.1.1]
    curl -s https://rest.kegg.jp/list/ko | \
        awk -F'\t' 'BEGIN{OFS="\t"} {gsub(/^ko:/, "", $1); print $1, $2}' > ko_to_name.tsv
fi

# -----------------------------------------------------------------
# 2. TaxID -> Scientific Name (NCBI Taxonomy dump から取得)
# -----------------------------------------------------------------
if [ ! -f taxid_to_name.tsv ]; then
    echo "Downloading NCBI Taxonomy names..."
    if [ ! -f taxdump.tar.gz ]; then
        wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    fi
    
    # names.dmp のみ展開
    tar -zxf taxdump.tar.gz names.dmp
    
    # "scientific name" の行だけを抽出し、taxid と name の 2 カラム TSV に変換
    # names.dmp のフォーマット: tax_id\t|\tname_txt\t|\tunique name\t|\tname class\t|
    awk -F'\t\\|\t' 'BEGIN{OFS="\t"} $4 ~ /^scientific name/ {
        gsub(/\t\|/, "", $4);
        print $1, $2
    }' names.dmp > taxid_to_name.tsv
    
    # 一時ファイルの削除
    rm -f names.dmp taxdump.tar.gz
fi

popd
popd
