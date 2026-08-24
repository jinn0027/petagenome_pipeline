#!/bin/bash
set -euo pipefail

### uniprot_to_{kegggenes,refseq,taxid,eggnog,biocyc,react}

# 2026/06/10 version
if [ ! -f idmapping.dat.gz ] ; then
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
fi

zcat idmapping.dat.gz | awk -F'\t' '
{
    ac = $1;
    type = $2;
    val = $3;

    if (type == "KEGG") {
        print ac "\t" val > "uniprot_to_kegggenes.tsv";
    }
    else if (type == "RefSeq") {
        print ac "\t" val > "uniprot_to_refseq.tsv";
    }
    else if (type == "NCBI_TaxID") {
        print ac "\t" val > "uniprot_to_taxid.tsv";
    }
    else if (type == "eggNOG") {
        print ac "\t" val > "uniprot_to_eggnog.tsv";
    }
    else if (type == "BioCyc") {
        print ac "\t" val > "uniprot_to_biocyc.tsv";
    }
    else if (type == "Reactome") {
        print ac "\t" val > "uniprot_to_react.tsv";
    }
    else if (type == "PlantReactome") {
        print ac "\t" val > "uniprot_to_plantreact.tsv";
    }
    else if (type == "TCDB") {
        print ac "\t" val > "uniprot_to_tcdb.tsv";
    }
    else if (type == "PHI-base") {
        print ac "\t" val > "uniprot_to_phibase.tsv";
    }
    else if (type == "GuidetoPHARMACOLOGY") {
        print ac "\t" val > "uniprot_to_pharmacology.tsv";
    }
}'

awk 'BEGIN { FS="\t"; OFS="\t" } {
    # 2列目のTC番号を "." で分割して、上位3つを取得
    n = split($2, a, ".")
    if (n >= 3) {
        family_id = a[1] "." a[2] "." a[3]
        print $1, family_id
    }
}' uniprot_to_tcdb.tsv > uniprot_to_tcdb_lough.tsv

### uniprot_to_go.tsv

if [ ! -f goa_uniprot_all.gaf.gz ] ; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
fi

zcat goa_uniprot_all.gaf.gz | awk -F'\t' '!/^!/ {print $2 "\t" $5}' | sort -u > uniprot_to_go.tsv

echo "抽出が完了しました。生成されたファイル一覧:"
ls -lh uniprot_to_*.tsv
