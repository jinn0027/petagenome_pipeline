#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義
params.memory = 32
params.threads = 8

// 2. タスク固有の上限値
def QIIME2_GG2_WGS_MAX_MEMORY = 64
def QIIME2_GG2_WGS_MAX_THREADS = 32

params.qiime2_gg2_wgs_memory = Math.min(params.memory as Integer, QIIME2_GG2_WGS_MAX_MEMORY)
params.qiime2_gg2_wgs_threads = Math.min(params.threads as Integer, QIIME2_GG2_WGS_MAX_THREADS)

// 3. rrnDBのコピー数探索モード ('right': 右側から細かい階層へ遡る, 'genus': 属レベル固定)
params.qiime2_gg2_rrndb_mode = 'right'

// 4. GXXXX (G[0-9]) ノードのみにバックボーンを絞るかどうか（デフォルト: true）
params.qiime2_gg2_only_gxxxx = true

// 5. 機能アノテーション関連パラメータ
if (!params.containsKey('annotation_table')) {
    params.annotation_table = "${params.petagenomeDir}/data/greengenes2/genome_annotations_table.tsv"
}
params.functional_annotations = "KO,MetaCyc"

// リファレンスファイルのパスのデフォルト設定
if (!params.containsKey('petagenomeDir') || !params.petagenomeDir) {
    error "Error: 'petagenomeDir' parameter is not specified. Please provide it via command line or config."
}

params.qiime2_gg2_wgs_backbone_fna = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.full-length.fna.qza"
params.qiime2_gg2_wgs_taxonomy = "${params.petagenomeDir}/data/greengenes2/2024.09.backbone.tax.qza"
params.qiime2_gg2_wgs_rrndb_stats = "${params.petagenomeDir}/data/rrnDB/rrnDB-5.10_pantaxa_stats_RDP.tsv.gz"

include { createNullParamsChannel; getParam; clusterOptions; processProfile; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

// ==========================================
// Greengenes2 による集計プロセス (WGS用)
// ==========================================
process qiime2_greengenes2_wgs {
    tag "${pair_id} (WGS)"
    container = "${params.petagenomeDir}/modules/qiime2/qiime2.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.qiime2_gg2_wgs_memory}"
    def threads = "${params.qiime2_gg2_wgs_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"
    
    input:
        tuple val(p), val(pair_id), path(sam_file)
        path backbone_fna
        path taxonomy
        path rrndb_stats
        val rrndb_mode
        path annotation_table
        val target_annots

    output:
        tuple val(pair_id), 
            path("${pair_id}/feature-table.tsv"), 
            path("${pair_id}/feature_counts.tsv"),
            path("${pair_id}/representatives.fasta"),
            path("${pair_id}/taxonomy.tsv"),
            path("${pair_id}/taxonomy_counts.tsv"),
            path("${pair_id}/*_functional_counts.tsv")

    script:
        """
        export PYTHONWARNINGS="ignore"
        export XDG_CONFIG_HOME=/tmp/qiime2_config
        export MPLCONFIGDIR=/tmp/matplotlib_config
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        export FONTCONFIG_PATH=/tmp/fontconfig

        echo "${processProfile(task)}" | tee prof.txt
        mkdir -p ${pair_id}

        echo "rrnDB search mode: ${rrndb_mode}"

        # 1. Bowtie2のSAMファイルからリファレンス配列ごとのカウントを集計
        python3 -c "
import sys
from collections import Counter

counts = Counter()
sam_path = '${sam_file}'

with open(sam_path, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        parts = line.strip().split('\\t')
        if len(parts) > 2:
            ref_id = parts[2]
            if ref_id != '*':
                counts[ref_id] += 1

sample_id = '${pair_id}'
with open('${pair_id}/feature-table.tsv', 'w') as out:
    out.write(f'# Constructed from biom file\\n#OTU ID\\t{sample_id}\\n')
    for ref_id, count in counts.items():
         out.write(f'{ref_id}\\t{count}\\n')
"

        # 2. 代表配列 (.qza) から代表配列 FASTA をエクスポートする
        qiime tools export \
            --input-path ${backbone_fna} \
            --output-path exported_ref

        if [ -f exported_ref/dna-sequences.fasta ]; then
            cp exported_ref/dna-sequences.fasta ${pair_id}/representatives.fasta
        elif ls exported_ref/*.fasta 1> /dev/null 2>&1; then
            cp exported_ref/*.fasta ${pair_id}/representatives.fasta
        else
            cp exported_ref/*.fna ${pair_id}/representatives.fasta
        fi

        # 3. タクソノミ情報をTSVに変換（安全なPython処理）
        qiime tools export \
            --input-path ${taxonomy} \
            --output-path exported_taxonomy

        python3 -c "
import os, glob, shutil

exported_dir = 'exported_taxonomy'
target_dest = '${pair_id}/taxonomy.tsv'

candidates = [
    os.path.join(exported_dir, 'taxonomy.tsv'),
    os.path.join(exported_dir, 'consensus_assignments.tsv')
]

found = False
for path in candidates:
    if os.path.exists(path):
        shutil.copy(path, target_dest)
        found = True
        break

if not found:
    all_tsvs = glob.glob(os.path.join(exported_dir, '**', '*.tsv'), recursive=True)
    if all_tsvs:
        shutil.copy(all_tsvs[0], target_dest)
        found = True

if not found:
    raise FileNotFoundError(f'Taxonomy tsv file not found in {exported_dir}')
"

        # 4. feature-table.tsv と taxonomy.tsv の結合・rrnDB補正 (タクソノミ集計)
        python3 ${params.petagenomeDir}/scripts/Python/parse_taxonomy.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${pair_id}/taxonomy_counts.tsv \
            --mode ${rrndb_mode}

        # 4.5. フィーチャー単位のカウントテーブルに対しても rrnDB補正・正規化・ソートを行い feature_counts.tsv を出力する
        python3 -c "
import pandas as pd
import zipfile
import gzip
import os

feat_path = '${pair_id}/feature-table.tsv'
tax_path = '${pair_id}/taxonomy.tsv'
rrndb_path = '${rrndb_stats}'
mode = '${rrndb_mode}'
output_path = '${pair_id}/feature_counts.tsv'

df_feat = pd.read_csv(feat_path, sep='\t', skiprows=1, index_col=0)
df_tax = pd.read_csv(tax_path, sep='\t', index_col=0)

sample_cols = df_feat.select_dtypes(include=['number']).columns
df_counts = df_feat[sample_cols].copy()

rrndb_df = None
if rrndb_path.endswith('.zip'):
    with zipfile.ZipFile(rrndb_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith('.tsv') or filename.endswith('.txt'):
                with z.open(filename) as f:
                    rrndb_df = pd.read_csv(f, sep='\t')
                break
elif rrndb_path.endswith('.gz'):
    with gzip.open(rrndb_path, 'rt') as f:
        rrndb_df = pd.read_csv(f, sep='\t')
else:
    rrndb_df = pd.read_csv(rrndb_path, sep='\t')

copy_num_dict = {}
if rrndb_df is not None:
    genus_col = next((c for c in ['genus', 'tax_name', 'name'] if c in rrndb_df.columns), rrndb_df.columns[0])
    mean_col = next((c for c in ['mean', 'copy_number'] if c in rrndb_df.columns), rrndb_df.columns[-1])
    for _, row in rrndb_df.iterrows():
        tax_name = str(row[genus_col]).strip()
        try:
            copy_num_dict[tax_name] = float(row[mean_col])
        except ValueError:
            pass

def get_rrn_copy(tax_str):
    if pd.isna(tax_str):
        return 1.0
    parts = [p.strip() for p in str(tax_str).split(';')]
    if mode == 'genus' and len(parts) >= 6:
        g = parts[5].replace('g__', '')
        if g in copy_num_dict:
            return copy_num_dict[g]
    for part in reversed(parts):
        clean_name = part.replace('g__', '').replace('s__', '').replace('f__', '')
        if clean_name in copy_num_dict:
            return copy_num_dict[clean_name]
    return 1.0

copy_nums = pd.Series(1.0, index=df_counts.index)
if 'Taxon' in df_tax.columns:
    tax_series = df_tax['Taxon']
elif df_tax.shape[1] > 0:
    tax_series = df_tax.iloc[:, 0]
else:
    tax_series = pd.Series(index=df_counts.index, dtype=str)

for idx in df_counts.index:
    if idx in tax_series.index:
        copy_nums[idx] = get_rrn_copy(tax_series[idx])

df_corrected = df_counts.div(copy_nums, axis=0)
df_norm = df_corrected.div(df_corrected.sum(axis=0), axis=1)

df_norm['Total'] = df_norm.sum(axis=1)
df_sorted = df_norm.sort_values(by='Total', ascending=False).drop(columns=['Total'])
df_sorted.to_csv(output_path, sep='\t')
"

        # 5. 機能アノテーション集計
        python3 ${params.petagenomeDir}/scripts/Python/parse_functional_profiles.py \
            ${pair_id}/feature-table.tsv \
            ${pair_id}/taxonomy.tsv \
            ${rrndb_stats} \
            ${annotation_table} \
            "${target_annots}" \
            "${pair_id}" \
            --mode ${rrndb_mode}
        """
}

// ==========================================
// サブワークフロー定義
// ==========================================
workflow QIIME2_GREENGENES2_WGS_SUB {
    take:
    p
    sam_ch
    backbone_fna
    taxonomy
    rrndb_stats
    rrndb_mode
    annotation_table
    target_annots

    main:
    out = qiime2_greengenes2_wgs(
        sam_ch,
        backbone_fna,
        taxonomy,
        rrndb_stats,
        rrndb_mode,
        annotation_table,
        target_annots
    )

    emit:
    out = out
}