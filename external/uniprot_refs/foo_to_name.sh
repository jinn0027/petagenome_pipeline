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
        raw_ko = $7;
        raw_name = $8;

        if (raw_ko != "") {
            # KO の加工 (| 以降を削除)
            split(raw_ko, ko_arr, "|");
            ko = ko_arr[1];

            # Name の加工 (";" および "," でも分割して、個別の名前にバラす)
            if (raw_name != "") {
                # まず ";" で分割
                n_names = split(raw_name, name_parts, ";");
                for (j=1; j<=n_names; j++) {
                    # さらにカンマ "," が含まれている場合もあるので分割する
                    n_sub = split(name_parts[j], sub_parts, ",");
                    for (k=1; k<=n_sub; k++) {
                        # 前後の空白を除去しつつ "|以降" を削る
                        sub(/^[ \t]+|[ \t]+$/, "", sub_parts[k]);
                        sub(/\|.*$/, "", sub_parts[k]);
                        trimmed = sub_parts[k];

                        if (trimmed != "") {
                            # KOごとにユニークな名前を登録
                            if (!((ko, trimmed) in registered)) {
                                registered[ko, trimmed] = 1;
                            
                                # リストとして蓄積
                                if (ko in ko_name_count) {
                                    count = ++ko_name_count[ko];
                                    ko_name_array[ko, count] = trimmed;
                                } else {
                                    ko_name_count[ko] = 1;
                                    ko_name_array[ko, 1] = trimmed;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    END {
        for (ko in ko_name_count) {
            line = "";
            for (i=1; i<=ko_name_count[ko]; i++) {
                if (i == 1) {
                    line = ko_name_array[ko, i];
                } else {
                    line = line ";" ko_name_array[ko, i];
                }
            }
            print ko "\t" line;
        }
    }' > ko_to_name.tsv

fi

