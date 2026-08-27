#!/bin/bash

### go_to_name.tsv

if [ ! -f go-basic.obo ] ; then
    wget https://purl.obolibrary.org/obo/go/go-basic.obo
fi

if [ ! -f go_to_name.tsv ] ; then
    awk '
        BEGIN { id=""; name=""; namespace="" }
        /^\[Term\]/ { id=""; name=""; namespace="" }
        /^id: / { id=$2 }
        /^name: / { 
            # "name: " 以降の文字列をすべて取得する
            name=$0; 
            sub(/^name: /, "", name);
        }
        /^namespace: / { 
            namespace=$2; 
            # タームのブロックが終わったら出力
            if (id != "" && name != "") {
                print id "\t" name "\t" namespace;
            }
        }' go-basic.obo > go_to_name.tsv
fi

### taxid_to_name.tsv

if [ ! -f taxdump.tar.gz ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
fi

if [ ! -f taxid_to_name.tsv ] ; then

    tar -zxf taxdump.tar.gz names.dmp
    
    # "scientific name" の行だけを抽出し、taxid と name の 2 カラム TSV に変換
    # names.dmp のフォーマット: tax_id\t|\tname_txt\t|\tunique name\t|\tname class\t|
    awk -F'\t\\|\t' 'BEGIN{OFS="\t"} $4 ~ /^scientific name/ {
    	gsub(/\t\|/, "", $4);
    	print $1, $2
	}' names.dmp > taxid_to_name.tsv
    
    rm -f names.dmp
fi

### react_to_name.tsv

# dl page:   https://reactome.org/download-data
if [ ! -f ReactomePathways.txt ] ; then
    wget https://download.reactome.org/97/ReactomePathways.txt
fi

if [ ! -f react_to_name.tsv ] ; then
    awk -F '\t' '{print($1 "\t" $2)}' ReactomePathways.txt > react_to_name.tsv
fi

### tcdb_to_name.tsv

if [ ! -f tcdb_families.tsv ] ; then
    curl -s "https://www.tcdb.org/cgi-bin/projectv/public/families.py" > tcdb_families.tsv
fi

if [ ! -f tcdb_to_name.tsv ] ; then
    awk -F '\t' '{print($1 "\t" $2)}' tcdb_families.tsv > tcdb_to_name.tsv
fi

### ko_to_name.tsv

if [ ! -f e7.og_info_kegg_go.tsv.gz ] ; then
    wget https://eggnogdb.org/public/eggnog7/e7.og_info_kegg_go.tsv.gz
fi

if [ ! -f ko_to_name.tsv ] ; then
    zcat e7.og_info_kegg_go.tsv.gz | awk -F'\t' '{
        # $6がタンパク質リスト、$7がKO、$8がKO名
        proteins = $6;
        ko = $7;
        name = $8;
    
        if (ko != "" && name != "") {
            print ko "\t" name > "ko_to_name.tsv";
   	}
    }'
fi

