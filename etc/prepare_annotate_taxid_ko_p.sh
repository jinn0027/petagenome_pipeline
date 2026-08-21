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

if [ ! -f uniprot_to_taxid.tsv  ] || [ ! -f uniprot_to_refseq.tsv ] || [ ! -f uniprot_to_gene.tsv ] || [ ! -f uniprot_to_go.tsv ]; then
    if [ ! -f idmapping_selected.tab.gz ] ; then
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    fi
    zcat idmapping_selected.tab.gz | awk -F'\t' ' { print $1 "\t" $3 > "uniprot_to_taxid.tsv"; print $1 "\t" $4 > "uniprot_to_refseq.tsv"; print $1 "\t" $5 > "uniprot_to_gene.tsv" ; print $1 "\t" $6 > "uniprot_to_go.tsv" ; } '
    for tsv in $(ls uniprot_to_*.tsv)
    do
	sed -i 's#; #;#g' ${tsv}
    done
fi

if [ ! -f refseq_to_ko.tsv ] || [ ! -f ko_to_name.tsv ] ; then
    if [ ! -f e7.og_info_kegg_go.tsv.gz ] ; then
	wget https://eggnogdb.org/public/eggnog7/e7.og_info_kegg_go.tsv.gz
    fi
    
    zcat e7.og_info_kegg_go.tsv.gz | awk -F'\t' '{
        # $6がタンパク質リスト、$7がKO、$8がKO名
    	proteins = $6;
    	ko = $7;
	name = $8;
    
        if (proteins != "" && ko != "") {
            n = split(proteins, p_arr, ",");
            for (i=1; i<=n; i++) {
                # 必要に応じてここで RefSeq などのIDフォーマットにフィルタリング
                print p_arr[i] "\t" ko > "refseq_to_ko.tsv1";
                print ko "\t" name > "ko_to_name.tsv1";
            }
        }
    }'

    # 2. ソート・ユニーク化と整形
    for tsv in refseq_to_ko.tsv ko_to_name.tsv; do
        if [ -f "${tsv}1" ]; then
            sort -u "${tsv}1" > "${tsv}"
            sed -i 's#; #;#g' "${tsv}"
            sed -i 's#, #,#g' "${tsv}"
        fi
    done

    # 3. 一時ファイルの削除
    rm -f *.tsv1
fi


if [ ! -f uniprot_to_ko.tsv ] ; then
    python ${MY_DIR}/etc/join_mapping.py -ab uniprot_to_refseq.tsv -bc refseq_to_ko.tsv -o uniprot_to_ko.tsv
fi

# -----------------------------------------------------------------
# TaxID -> Scientific Name (NCBI Taxonomy dump から取得)
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
