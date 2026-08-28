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
            split(raw_ko, ko_arr, "|");
            ko = ko_arr[1];
    
            if (raw_name != "") {
                n_names = split(raw_name, name_parts, ";");
                for (j=1; j<=n_names; j++) {
                    n_sub = split(name_parts[j], sub_parts, ",");
                    for (k=1; k<=n_sub; k++) {
                        sub(/^[ \t]+|[ \t]+$/, "", sub_parts[k]);
                        sub(/\|.*$/, "", sub_parts[k]);
                        trimmed = sub_parts[k];
    
                        if (trimmed != "") {
                            # 出現回数をカウント
                            freq[ko, trimmed]++;
                            
                            # 初めて登場した (ko, trimmed) の組み合わせであれば、リスト用の配列に蓄積
                            if (!((ko, trimmed) in seen)) {
                                seen[ko, trimmed] = 1;
                                
                                # カウントを保持して後で配列化できるようにする
                                idx = ++ko_count[ko];
                                ko_elements[ko, idx] = trimmed;
                            }
                        }
                    }
                }
            }
        }
    }
    END {
        for (ko in ko_count) {
            n = ko_count[ko];
            
            # 配列から一時的な names 配列に移してソート
            for (i = 1; i <= n; i++) {
                names[i] = ko_elements[ko, i];
            }
    
            # 出現回数（freq）の降順でソート（同率の場合は名前順）
            for (i = 1; i <= n; i++) {
                for (j = i + 1; j <= n; j++) {
                    f_i = freq[ko, names[i]];
                    f_j = freq[ko, names[j]];
                    if (f_j > f_i || (f_j == f_i && names[j] < names[i])) {
                        temp = names[i];
                        names[i] = names[j];
                        names[j] = temp;
                    }
                }
            }
    
            # 頻度に応じて ">" または "=" で結合
            line = names[1];
            for (i = 2; i <= n; i++) {
                prev_freq = freq[ko, names[i-1]];
                curr_freq = freq[ko, names[i]];
                
                if (curr_freq == prev_freq) {
                    line = line "=" names[i];
                } else {
                    line = line ">" names[i];
                }
            }
            print ko "\t" line;
        }
    }' > ko_to_name.tsv
fi

