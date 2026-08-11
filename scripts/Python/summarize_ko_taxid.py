#!/usr/bin/env python3
import sys
import os
import argparse
import pandas as pd

def load_name_map(map_filepath, id_col_name):
    """
    ヘッダーなしTSV (1列目: ID, 2列目: Name) を読み込んで DataFrame として返す。
    ファイルが存在しない場合や None の場合は空の DataFrame を返す。
    """
    if not map_filepath or not os.path.exists(map_filepath):
        return None
    
    try:
        df_map = pd.read_csv(
            map_filepath, 
            sep='\t', 
            header=None, 
            names=[id_col_name, 'name'], 
            dtype=str
        )
        # 重複IDがある場合は最初の定義を採用
        df_map = df_map.drop_duplicates(subset=[id_col_name])
        return df_map
    except Exception as e:
        sys.stderr.write(f"Warning: Failed to load name map file {map_filepath}: {e}\n")
        return None

def process_sample(filepath, out_prefix, min_anno_pident=0.0, ko_map_path=None, taxid_map_path=None):
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # 低相同性フィルタリング
    if min_anno_pident > 0.0 and 'anno_pident' in df.columns:
        df['anno_pident_num'] = pd.to_numeric(df['anno_pident'], errors='coerce')
        df = df[df['anno_pident_num'] >= min_anno_pident]

    # マッピングファイルのロード
    ko_map_df = load_name_map(ko_map_path, 'ko')
    taxid_map_df = load_name_map(taxid_map_path, 'taxid')

    # --- 1. KO 集計 ---
    if 'ko' in df.columns:
        df_ko = df[['orf_id', 'ko']].dropna().drop_duplicates()
        df_ko = df_ko[df_ko['ko'] != 'N/A']
        
        total_ko_orfs = len(df_ko['orf_id'].unique())
        ko_counts = df_ko['ko'].value_counts().reset_index()
        ko_counts.columns = ['ko', 'count']
        ko_counts['ratio'] = ko_counts['count'] / total_ko_orfs if total_ko_orfs > 0 else 0.0

        # 名称のマージ
        if ko_map_df is not None:
            ko_counts = pd.merge(ko_counts, ko_map_df, on='ko', how='left')
            ko_counts['name'] = ko_counts['name'].fillna('Unknown')
            # 列順の整理: ko, name, count, ratio
            ko_counts = ko_counts[['ko', 'name', 'count', 'ratio']]
    else:
        ko_counts = pd.DataFrame(columns=['ko', 'count', 'ratio'])

    # --- 2. TaxID 集計 ---
    if 'taxid' in df.columns:
        df_tax = df[['orf_id', 'taxid']].dropna().drop_duplicates()
        df_tax = df_tax[df_tax['taxid'] != 'N/A'].copy()
        df_tax['taxid'] = df_tax['taxid'].str.split(';')
        df_tax = df_tax.explode('taxid')
        df_tax['taxid'] = df_tax['taxid'].str.strip()
        df_tax = df_tax.drop_duplicates()
        
        total_tax_orfs = len(df_tax['orf_id'].unique())
        tax_counts = df_tax['taxid'].value_counts().reset_index()
        tax_counts.columns = ['taxid', 'count']
        tax_counts['ratio'] = tax_counts['count'] / total_tax_orfs if total_tax_orfs > 0 else 0.0

        # 名称のマージ
        if taxid_map_df is not None:
            tax_counts = pd.merge(tax_counts, taxid_map_df, on='taxid', how='left')
            tax_counts['name'] = tax_counts['name'].fillna('Unknown')
            # 列順の整理: taxid, name, count, ratio
            tax_counts = tax_counts[['taxid', 'name', 'count', 'ratio']]
    else:
        tax_counts = pd.DataFrame(columns=['taxid', 'count', 'ratio'])

    # 保存
    ko_counts.to_csv(f"{out_prefix}.ko_summary.tsv", sep='\t', index=False)
    tax_counts.to_csv(f"{out_prefix}.taxid_summary.tsv", sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize KO and TaxID ratios per sample.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument("-p", "--min_pident", type=float, default=0.0, help="Minimum anno_pident cutoff")
    parser.add_argument("--ko_name_map", default=None, help="Path to ko_to_name.tsv (headerless: KO, Name)")
    parser.add_argument("--taxid_name_map", default=None, help="Path to taxid_to_name.tsv (headerless: TaxID, Name)")

    args = parser.parse_args()
    process_sample(
        args.input, 
        args.out_prefix, 
        args.min_pident, 
        ko_map_path=args.ko_name_map, 
        taxid_map_path=args.taxid_name_map
    )