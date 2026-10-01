#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import gzip
import zipfile
import os

def open_table(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    elif path.endswith(".zip"):
        z = zipfile.ZipFile(path)
        namelist = z.namelist()
        return z.open(namelist[0], "r")
    else:
        return open(path, "r")

def load_rrndb(rrndb_stats_path, mode="right"):
    rrndb_dict = {}
    with open_table(rrndb_stats_path) as f:
        for line in f:
            if isinstance(line, bytes):
                line = line.decode('utf-8', errors='ignore')
                
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                tax_str = parts[0]
                try:
                    mean_rrn = float(parts[1])
                    rrndb_dict[tax_str] = mean_rrn
                except ValueError:
                    continue
    return rrndb_dict

def get_rrn_copy_number(taxonomy_str, rrndb_dict, mode="right"):
    if not isinstance(taxonomy_str, str):
        return 1.0

    levels = [l.strip() for l in taxonomy_str.split(";")]
    cleaned_levels = []
    for l in levels:
        if "__" in l:
            cleaned_levels.append(l.split("__", 1)[1])
        else:
            cleaned_levels.append(l)

    if mode == "genus":
        for i in range(min(len(cleaned_levels), 6) - 1, -1, -1):
            sub_tax = ";".join(cleaned_levels[:i+1])
            if sub_tax in rrndb_dict:
                return rrndb_dict[sub_tax]
    else:  # 'right' モード: 最も細かい階層から順に遡る
        for i in range(len(cleaned_levels), 0, -1):
            sub_tax = ";".join(cleaned_levels[:i])
            if sub_tax in rrndb_dict:
                return rrndb_dict[sub_tax]

    return 1.0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("feature_table", help="Path to feature-table.tsv")
    parser.add_argument("taxonomy_file", help="Path to taxonomy.tsv")
    parser.add_argument("rrndb_stats", help="Path to rrnDB stats")
    parser.add_argument("annotation_table", help="Path to genome_annotations_table.tsv")
    parser.add_argument("target_annots", help="Comma-separated annotations (e.g. KO,MetaCyc)")
    parser.add_argument("output_prefix", help="Output prefix directory/name")
    parser.add_argument("--mode", default="right", help="rrnDB mode")
    args = parser.parse_args()

    target_cols = [c.strip() for c in args.target_annots.split(",")]

    # 1. アノテーションテーブルの読み込み
    ann_df = pd.read_csv(args.annotation_table, sep="\t")
    if "GenomeID" in ann_df.columns:
        ann_df.set_index("GenomeID", inplace=True)
    elif "FeatureID" in ann_df.columns:
        ann_df.set_index("FeatureID", inplace=True)

    # 2. フィーチャーテーブルとタクソノミの読み込み（1サンプル列を前提に特定）
    ft_df = pd.read_csv(args.feature_table, sep="\t", skiprows=1)
    ft_df.rename(columns={ft_df.columns[0]: "FeatureID"}, inplace=True)
    
    # FeatureID以外の最初の列をサンプル名として直接取得
    sample_col = [c for c in ft_df.columns if c != "FeatureID"][0]
    
    tax_df = pd.read_csv(args.taxonomy_file, sep="\t")
    tax_id_col = tax_df.columns[0]
    tax_col = [c for c in tax_df.columns if 'tax' in c.lower() or 'consensus' in c.lower()][0]
    tax_df = tax_df[[tax_id_col, tax_col]]
    tax_df.columns = ["FeatureID", "Taxonomy"]

    # 3. rrnDBの読み込みとコピー数補正
    rrndb_dict = load_rrndb(args.rrndb_stats, mode=args.mode)
    merged_df = pd.merge(ft_df, tax_df, on="FeatureID", how="inner")

    # rrnDBコピー数による割り算
    corrected_data = []
    for idx, row in merged_df.iterrows():
        tax_str = row["Taxonomy"]
        copy_num = get_rrn_copy_number(tax_str, rrndb_dict, mode=args.mode)
        
        new_row = row.copy()
        new_row[sample_col] = float(row[sample_col]) / copy_num
        corrected_data.append(new_row)

    corr_df = pd.DataFrame(corrected_data)

    # 全体和正規化（サンプルの合計が1になるようにする）
    col_sum = corr_df[sample_col].sum()
    if col_sum > 0:
        corr_df[sample_col] = corr_df[sample_col] / col_sum

    # 4. アノテーションテーブルとの結合
    if corr_df["FeatureID"].isin(ann_df.index).any():
        corr_df = corr_df.join(ann_df, on="FeatureID", how="inner")
    else:
        corr_df = pd.merge(corr_df, ann_df, left_on="FeatureID", right_index=True, how="inner")

    os.makedirs(args.output_prefix, exist_ok=True)

    # 5. 各ターゲットアノテーションごとに展開・集計・ソート・出力
    for annot in target_cols:
        output_file = os.path.join(args.output_prefix, f"{annot.lower()}_functional_counts.tsv")
        
        if annot not in corr_df.columns:
            print(f"Warning: Annotation column '{annot}' not found in merged table.", file=sys.stderr)
            with open(output_file, "w") as f:
                f.write(f"FeatureID\t{sample_col}\n")
            continue

        # カンマ区切りの文字列を展開
        df_expanded = corr_df.copy()
        df_expanded[annot] = df_expanded[annot].astype(str).str.split(",")
        df_expanded = df_expanded.explode(annot)
        df_expanded[annot] = df_expanded[annot].str.strip()
        
        # 無効値の除外
        df_expanded = df_expanded[df_expanded[annot].notna() & (df_expanded[annot] != "") & (df_expanded[annot] != "nan")]

        # 機能コードごとにグループ集計して降順ソート
        grouped_df = df_expanded.groupby(annot)[sample_col].sum()
        sorted_grouped_df = grouped_df.sort_values(ascending=False)

        sorted_grouped_df.to_csv(output_file, sep="\t", header=True, index_label=annot)

    print(f"Functional profiles for {target_cols} calculated and saved successfully.")

if __name__ == "__main__":
    main()