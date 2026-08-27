#!/bin/bash
set -euo pipefail

### uniprot_to_{kegggenes,refseq,taxid,eggnog,biocyc,react}

# 2026/06/10 version
if [ ! -f idmapping.dat.gz ] ; then
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
fi

if [ ! -f uniprot_to_tcdb.tsv ] ; then
    
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

fi

if [ ! -f uniprot_to_tcdb_lough.tsv ] ; then
    awk 'BEGIN { FS="\t"; OFS="\t" } {
        # 2列目のTC番号を "." で分割して、上位3つを取得
        n = split($2, a, ".")
        if (n >= 3) {
            family_id = a[1] "." a[2] "." a[3]
            print $1, family_id
        }
    }' uniprot_to_tcdb.tsv > uniprot_to_tcdb_lough.tsv
fi

### uniprot_to_go.tsv

if [ ! -f goa_uniprot_all.gaf.gz ] ; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
fi

if [ ! -f uniprot_to_go.tsv ] ; then
    zcat goa_uniprot_all.gaf.gz | awk -F'\t' '!/^!/ {print $2 "\t" $5}' | sort -u > uniprot_to_go.tsv
fi

### uniprot_to_ko.tsv

if [ ! -f uniprot_ids.txt ] ; then
    zcat idmapping.dat.gz | awk '{print($1)}' | uniq > uniprot_ids.txt
fi

if [ ! -f uniprot_to_ko.tsv ] ; then
python <<EOF
import gzip
import sys

# 1. 1.5億行の uniprot_ids.txt をメモリ上のハッシュセットにロード
print("Loading uniprot_ids.txt...", file=sys.stderr)
uniprot_set = set()
with open("./uniprot_ids.txt", "r") as f:
  for line in f:
    uid = line.strip()
    if uid:
      uniprot_set.add(uid)
print(f"Loaded {len(uniprot_set)} IDs.", file=sys.stderr)

# 2. e7.og_info_kegg_go.tsv.gz を読み込みながら突合して出力
print("Processing e7.og_info_kegg_go.tsv.gz...", file=sys.stderr)
output_file = "uniprot_to_ko.tsv"

with gzip.open("e7.og_info_kegg_go.tsv.gz", "rt") as fin, open(
        output_file, "w"
) as fout:
  for line in fin:
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
      continue

    proteins = cols[5]  # $6
    raw_ko = cols[6]  # $7

    if proteins and raw_ko:
      ko = raw_ko.split("|")[0]
      for p in proteins.split(","):
        sub_parts = p.split(".")
        if len(sub_parts) == 2:
          prot_id = sub_parts[1]
          # O(1) で1.5億行のセットから判定
          if prot_id in uniprot_set:
            fout.write(f"{prot_id}\t{ko}\n")

print("Done!", file=sys.stderr)
EOF

fi

echo "抽出が完了しました。生成されたファイル一覧:"
ls -lh uniprot_to_*.tsv

