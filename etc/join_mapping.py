import argparse
import sys
from collections import defaultdict

def parse_and_add(val_str, target_set):
    """セミコロンやカンマ区切りの文字列をバラしてセットに追加する"""
    if not val_str or val_str == "-":
        return
    for item in val_str.replace(",", ";").split(";"):
        item = item.strip()
        if item:
            target_set.add(item)

def build_mapping(file_path):
    """A -> [Bの集合] の辞書を作る（複数行・セミコロン区切りに対応）"""
    mapping = defaultdict(set)
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            key = parts[0].strip()
            val_field = parts[1].strip()
            
            items = set()
            parse_and_add(val_field, items)
            
            for item in items:
                mapping[key].add(item)
    return mapping

def main():
    parser = argparse.ArgumentParser(description="Join two mapping TSV files (A->B and B->C) into A->C.")
    parser.add_argument("-ab", "--a_to_b", required=True, help="Path to A_to_B TSV file")
    parser.add_argument("-bc", "--b_to_c", required=True, help="Path to B_to_C TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to output A_to_C TSV file")
    
    args = parser.parse_args()
    
    print(f"Loading A -> B mapping from {args.a_to_b}...")
    a_to_b = build_mapping(args.a_to_b)
    
    print(f"Loading B -> C mapping from {args.b_to_c}...")
    b_to_c = build_mapping(args.b_to_c)
    
    print("Computing A -> C mapping...")
    a_to_c = defaultdict(set)
    
    for a, b_set in a_to_b.items():
        for b in b_set:
            # スコア付き（例: K00001|100.00）を考慮してクレンジング
            b_clean = b.split('|')[0] if '|' in b else b
            
            if b_clean in b_to_c:
                for c in b_to_c[b_clean]:
                    a_to_c[a].add(c)
            if b in b_to_c:
                for c in b_to_c[b]:
                    a_to_c[a].add(c)

    print(f"Writing to {args.output}...")
    with open(args.output, 'w', encoding='utf-8') as f:
        for a in sorted(a_to_c.keys()):
            c_concatenated = ";".join(sorted(a_to_c[a]))
            f.write(f"{a}\t{c_concatenated}\n")
            
    print("Done!")

if __name__ == "__main__":
    main()
