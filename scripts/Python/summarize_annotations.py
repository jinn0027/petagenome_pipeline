#!/usr/bin/env python3
import sys
import os
import argparse
import pandas as pd

def load_name_map(map_filepath, id_col_name):
    """
    ヘッダーなしTSV (1列目: ID, 2列目: Name) を読み込んで DataFrame として返す。
    ファイルが存在しない場合や None の場合は None を返す。
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

def process_sample(filepath, out_prefix, min_anno_pident=0.0, target_attrs=None, name_maps=None):
    if target_attrs is None:
        target_attrs = ['taxid', 'ko']
    if name_maps is None:
        name_maps = {}

    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # 低相同性フィルタリング
    if min_anno_pident > 0.0 and 'anno_pident' in df.columns:
        df['anno_pident_num'] = pd.to_numeric(df['anno_pident'], errors='coerce')
        df = df[df['anno_pident_num'] >= min_anno_pident]

    # 各属性ごとにループして集計
    for attr in target_attrs:
        out_filename = f"{out_prefix}.{attr}.summary.tsv"
        
        if attr in df.columns:
            df_attr = df[['orf_id', attr]].dropna().drop_duplicates()
            df_attr = df_attr[df_attr[attr] != 'N/A'].copy()
            
            # セミコロン区切りの値を分割して展開
            df_attr[attr] = df_attr[attr].str.split(';')
            df_attr = df_attr.explode(attr)
            df_attr[attr] = df_attr[attr].str.strip()
            df_attr = df_attr.drop_duplicates()
            
            total_orfs = len(df_attr['orf_id'].unique())
            counts = df_attr[attr].value_counts().reset_index()
            counts.columns = [attr, 'count']
            counts['ratio'] = counts['count'] / total_orfs if total_orfs > 0 else 0.0

            # 該当属性の名前マップパスを取得してロード
            map_path = name_maps.get(attr)
            map_df = load_name_map(map_path, attr)

            if map_df is not None:
                # マップが存在する場合はマージして name 列を追加
                counts = pd.merge(counts, map_df, on=attr, how='left')
                counts['name'] = counts['name'].fillna('Unknown')
                # 列順の整理: [attr, 'name', 'count', 'ratio']
                counts = counts[[attr, 'name', 'count', 'ratio']]
            else:
                # マップが存在しない（未指定の）場合は name 列を含めない
                counts = counts[[attr, 'count', 'ratio']]
        else:
            # カラム自体が存在しない場合は空の DataFrame を作成
            map_path = name_maps.get(attr)
            map_df = load_name_map(map_path, attr)
            if map_df is not None:
                counts = pd.DataFrame(columns=[attr, 'name', 'count', 'ratio'])
            else:
                counts = pd.DataFrame(columns=[attr, 'count', 'ratio'])

        # 保存
        counts.to_csv(out_filename, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize annotation ratios per sample dynamically.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument("-p", "--min_pident", type=float, default=0.0, help="Minimum anno_pident cutoff")
    parser.add_argument("--target_attrs", required=True, help="Comma-separated list of target attributes (e.g., taxid,ko)")

    # 既知の属性だけでなく、将来的な拡張を見据えて動的に引数を受け付けられるようにする
    # Nextflow 側からは `--<attr>_name_map` という形式で渡される
    known_args, unknown = parser.parse_known_args()
    
    target_attrs = [attr.strip() for attr in known_args.target_attrs.split(',') if attr.strip()]
    
    # 動的に渡されたオプション引数から name_map を抽出する
    # 例: --ko_name_map path/to/ko.tsv -> name_maps['ko'] = 'path/to/ko.tsv'
    name_maps = {}
    i = 0
    while i < len(unknown):
        arg = unknown[i]
        if arg.startswith('--') and '_name_map' in arg:
            attr_name = arg.replace('--', '').replace('_name_map', '')
            if i + 1 < len(unknown):
                name_maps[attr_name] = unknown[i + 1]
                i += 2
                continue
        i += 1

    process_sample(
        known_args.input, 
        known_args.out_prefix, 
        known_args.min_pident, 
        target_attrs=target_attrs,
        name_maps=name_maps
    )