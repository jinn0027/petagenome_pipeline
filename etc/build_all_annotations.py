import lzma
import json
import csv
from collections import defaultdict

def build_all_genome_annotations(coords_path, uniref_map_path, map_paths_dict):
    """
    coords_path: proteins/coords.txt.xz
    uniref_map_path: function/uniref/uniref.map.xz
    map_paths_dict: 取得したい機能の名称とファイルパスの辞書
    """
    print("1. Loading protein coordinates (Genome -> Protein ID)...")
    protein_to_genome = {}
    current_genome = None
    
    with lzma.open(coords_path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_genome = line[1:].strip()
            elif current_genome:
                parts = line.split("\t")
                protein_id = f"{current_genome}_{parts[0]}"
                protein_to_genome[protein_id] = current_genome

    print(f"   Loaded {len(protein_to_genome):,} proteins.")

    print("2. Loading UniRef mapping (Protein ID -> UniRef)...")
    protein_to_uniref = {}
    with lzma.open(uniref_map_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protein_to_uniref[parts[0]] = parts[1]

    print(f"   Loaded {len(protein_to_uniref):,} UniRef mappings.")

    # 3. 各種機能マップのロード
    loaded_maps = {}
    for name, path in map_paths_dict.items():
        print(f"3. Loading {name} mapping from {path}...")
        feature_map = {}
        with lzma.open(path, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    key_id = parts[0]
                    vals = parts[1:]
                    feature_map[key_id] = vals
        loaded_maps[name] = feature_map
        print(f"   Loaded {len(feature_map):,} entries for {name}.")

    print("4. Chaining mappings to build genome-level annotations for all categories...")
    genome_annotations = defaultdict(lambda: defaultdict(set))
    
    matched_counts = {name: 0 for name in map_paths_dict}

    for prot_id, genome_id in protein_to_genome.items():
        uniref_id = protein_to_uniref.get(prot_id)
        
        for name, f_map in loaded_maps.items():
            keys_to_check = [uniref_id, prot_id]
            matched = False
            for k in keys_to_check:
                if k and k in f_map:
                    for val in f_map[k]:
                        genome_annotations[genome_id][name].add(val)
                    matched = True
                    break
            if matched:
                matched_counts[name] += 1

    print("\n   Matching summary:")
    for name, count in matched_counts.items():
        print(f"     - {name}: {count:,} proteins matched")

    result = {}
    for g_id, ann_dict in genome_annotations.items():
        result[g_id] = {name: sorted(list(vals)) for name, vals in ann_dict.items()}

    print(f"\n   Done! Total genomes with annotations: {len(result):,}")
    return result

def save_annotations_to_tsv(all_annotations, output_path):
    """
    ゲodムごとに、各種アノテーションをタブ区切り（リストはカンマ結合）でファイル出力する
    """
    print(f"Saving results to TSV: {output_path}...")
    
    # どのアノテーションカテゴリが存在するか動的に収集
    categories = set()
    for g_data in all_annotations.values():
        categories.update(g_data.keys())
    categories = sorted(list(categories))
    
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        # ヘッダー行
        header = ["GenomeID"] + categories
        writer.writerow(header)
        
        # 各行（ゲノムごと）
        for genome_id, data in all_annotations.items():
            row = [genome_id]
            for cat in categories:
                # リストをカンマ区切りの文字列にして格納（空の場合は空文字）
                items = data.get(cat, [])
                row.append(",".join(items))
            writer.writerow(row)
    print("   TSV export completed.")

def save_annotations_to_json(all_annotations, output_path):
    """
    Pythonの辞書構造のままJSON形式で保存する（後からの再利用に便利）
    """
    print(f"Saving results to JSON: {output_path}...")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(all_annotations, f, ensure_ascii=False, indent=2)
    print("   JSON export completed.")

if __name__ == "__main__":
    coords_file = "./proteins/coords.txt.xz"
    uniref_file = "./function/uniref/uniref.map.xz"

    target_maps = {
        "KO": "./function/kegg/ko.map.xz",
        "MetaCyc": "./function/metacyc/protein.map.xz",
        "GO_All": "./function/go/all.map.xz",
        "GO_Process": "./function/go/process.map.xz",
        "eggNOG": "./function/eggnog/eggnog.map.xz",
        "OrthoDB": "./function/orthodb/orthodb.map.xz",
        "RefSeq": "./function/refseq/refseq.map.xz"
    }

    # アノテーションの構築
    all_annotations = build_all_genome_annotations(coords_file, uniref_file, target_maps)

    if all_annotations:
        # 1. TSV形式で保存（エクセルや他のツール、R・Pythonでの読み込み用）
        save_annotations_to_tsv(all_annotations, "genome_annotations_table.tsv")

        # 2. JSON形式で保存（プログラムで構造を保持したまま読み込む用）
        #save_annotations_to_json(all_annotations, "genome_annotations_dict.json")

        # サンプルの表示
        sample_genome = list(all_annotations.keys())[0]
        print(f"\n[Sample] Genome ID: {sample_genome}")
        for ann_type, items in all_annotations[sample_genome].items():
            print(f"  - {ann_type} (total {len(items)}): {items[:5]} ...")