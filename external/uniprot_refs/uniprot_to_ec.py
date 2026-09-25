import re
import argparse

def parse_enzyme_dat_to_map(enzyme_dat_path, out_path):
    # UniProt ID ごとに複数の EC 番号を格納するための辞書
    ac_to_ecs = {}
    
    current_ec = None
    current_acs = []
    
    print(f"Reading {enzyme_dat_path}...")
    with open(enzyme_dat_path, 'r', encoding='utf-8') as fin:
        for line in fin:
            line = line.rstrip('\r\n')
            
            if line.startswith('ID   '):
                current_ec = line[5:].strip()
                current_acs = []
            elif line.startswith('DR   '):
                # DR行からUniProtのアクセッションを抽出
                parts = line[5:].split(';')
                for part in parts:
                    tokens = part.strip().split(',')
                    if tokens:
                        ac = tokens[0].strip()
                        # 簡易的なUniProtアクセッションの形式チェック
                        if re.match(r'^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-V][0-9][A-Z0-9]{2}[0-9]([A-Z0-9][A-Z0-9]{2}[0-9])?$', ac):
                            current_acs.append(ac)
            elif line == '//':
                if current_ec and current_acs:
                    for ac in current_acs:
                        if ac not in ac_to_ecs:
                            ac_to_ecs[ac] = set()
                        ac_to_ecs[ac].add(current_ec)
                current_ec = None
                current_acs = []
                
    # ファイルに出力（複数ある場合はセミコロン区切り）
    count = 0
    with open(out_path, 'w', encoding='utf-8') as fout:
        for ac, ecs in sorted(ac_to_ecs.items()):
            ec_str = ";".join(sorted(ecs))
            fout.write(f"{ac}\t{ec_str}\n")
            count += 1
            
    print(f"Done! Extracted {count} UniProt-to-EC mappings to {out_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract UniProt-to-EC mapping from enzyme.dat with semicolon separation.")
    parser.add_argument("-i", "--input", default="enzyme.dat", help="Input enzyme.dat file")
    parser.add_argument("-o", "--output", default="uniprot_to_ec.tsv", help="Output TSV file")
    
    args = parser.parse_args()
    parse_enzyme_dat_to_map(args.input, args.output)
