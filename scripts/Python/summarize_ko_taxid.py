#!/usr/bin/env python3
import sys
import argparse
import pandas as pd

def process_sample(filepath, out_prefix, min_anno_pident=0.0):
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # 低相同性フィルタリング
    if min_anno_pident > 0.0 and 'anno_pident' in df.columns:
        df['anno_pident_num'] = pd.to_numeric(df['anno_pident'], errors='coerce')
        df = df[df['anno_pident_num'] >= min_anno_pident]

    # --- 1. KO 集計 ---
    if 'ko' in df.columns:
        df_ko = df[['orf_id', 'ko']].dropna().drop_duplicates()
        df_ko = df_ko[df_ko['ko'] != 'N/A']
        
        total_ko_orfs = len(df_ko['orf_id'].unique())
        ko_counts = df_ko['ko'].value_counts().reset_index()
        ko_counts.columns = ['ko', 'count']
        ko_counts['ratio'] = ko_counts['count'] / total_ko_orfs if total_ko_orfs > 0 else 0.0
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

    args = parser.parse_args()
    process_sample(args.input, args.out_prefix, args.min_pident)
