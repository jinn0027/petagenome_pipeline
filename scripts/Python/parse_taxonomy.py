#!/usr/bin/env python3
import sys
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description="Parse taxonomy and normalize by rrnDB copy number.")
parser.add_argument("table_path", help="Path to feature-table.tsv")
parser.add_argument("taxonomy_path", help="Path to taxonomy.tsv")
parser.add_argument("rrndb_path", help="Path to rrnDB stats file")
parser.add_argument("output_path", help="Path to output taxonomy_counts.tsv")
parser.add_argument("--mode", default="right", choices=["right", "genus"], help="rrnDB copy number search mode ('right': fallback from fine to coarse, 'genus': genus level only)")

args = parser.parse_args()

table_path = args.table_path
taxonomy_path = args.taxonomy_path
rrndb_path = args.rrndb_path
output_path = args.output_path
rrndb_mode = args.mode

print(f"Generating copy-number normalized and sorted taxonomy_counts.tsv (mode: {rrndb_mode})...")

# rrnDB 統計ファイルの読み込み
rrndb_df = pd.read_csv(rrndb_path, sep="\t", compression="infer")
copy_dict = dict(zip(rrndb_df["name"], rrndb_df["mean"]))

matched_count = 0
unmatched_count = 0
unmatched_taxa_examples = []

def clean_taxon_name(name):
    # 余分なプレフィックス (g__, f__, o__ など) を削除
    name = re.sub(r"^[a-z]__", "", name).strip()
    # Greengenes2などで見られる末尾の付加コードなどをクレンジング
    name = re.sub(r"_[A-Z]_[0-9]+", "", name)
    name = re.sub(r"_[0-9]+$", "", name)
    return name

def get_copy_number(tax_str):
    global matched_count, unmatched_count, unmatched_taxa_examples
    
    # 1. 属レベル固定モード ('genus')
    if rrndb_mode == "genus":
        match_genus = re.search(r"g__([^;]+)", tax_str)
        if match_genus:
            g_name = match_genus.group(1)
            g_clean = clean_taxon_name(g_name)
            if g_clean in copy_dict and g_clean:
                matched_count += 1
                return copy_dict[g_clean]
    
    # 2. 右側（細かい階層）から遡るモード ('right')
    elif rrndb_mode == "right":
        # セミコロンで分割し、有効な階層を後ろから順に評価
        parts = [p.strip() for p in tax_str.split(";") if p.strip()]
        for part in reversed(parts):
            cleaned = clean_taxon_name(part)
            if not cleaned or cleaned.lower() in ["unassigned", "incertae_sedis", "uncultured"]:
                continue
            if cleaned in copy_dict:
                matched_count += 1
                return copy_dict[cleaned]
    
    # マッチしなかった場合
    unmatched_count += 1
    if len(unmatched_taxa_examples) < 10:
        unmatched_taxa_examples.append(tax_str)
    return 1.0

# feature-table.tsv のスキップ行特定と読み込み
skiprows = 0
with open(table_path, "r") as f:
    for i, line in enumerate(f):
        if line.startswith("#OTU ID") or line.startswith("#Feature ID"):
            skiprows = i
            break

table_df = pd.read_csv(table_path, sep="\t", skiprows=skiprows, index_col=0)

if table_df.index.name and table_df.index.name.startswith("#"):
    table_df.index.name = table_df.index.name.lstrip("#").strip()

tax_df = pd.read_csv(taxonomy_path, sep="\t", index_col=0)
tax_col = tax_df.columns[0]
tax_dict = tax_df[tax_col].to_dict()

counts = table_df.sum(axis=1)

tax_counts = {}
for feat_id, count in counts.items():
    taxon = tax_dict.get(str(feat_id), "k__Unassigned; p__; c__; o__; f__; g__; s__")
    tax_counts[taxon] = tax_counts.get(taxon, 0) + count

sample_name = table_df.columns[0] if len(table_df.columns) > 0 else "count"

summary_out = pd.DataFrame({
    "Taxonomy": list(tax_counts.keys()),
    sample_name: list(tax_counts.values())
})

summary_out["copy_number"] = summary_out["Taxonomy"].apply(get_copy_number)
summary_out["norm_val"] = summary_out[sample_name] / summary_out["copy_number"]

print("--- rrnDB Copy Number Matching Debug ---")
print(f"  Matched (found in rrnDB): {matched_count}")
print(f"  Unmatched (defaulted to 1.0): {unmatched_count}")
if unmatched_taxa_examples:
    print("  Unmatched taxa examples (up to 10):")
    for ex in unmatched_taxa_examples:
        print(f"    - {ex}")
print("----------------------------------------")

total_norm = summary_out["norm_val"].sum()
if total_norm > 0:
    summary_out[sample_name] = summary_out["norm_val"] / total_norm
else:
    summary_out[sample_name] = 0.0

summary_out = summary_out[["Taxonomy", sample_name]].sort_values(by=sample_name, ascending=False)
summary_out.to_csv(output_path, sep="\t", index=False)
print("taxonomy_counts.tsv generated, normalized, and sorted successfully.")