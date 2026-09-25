#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. 全体デフォルト値の定義（未定義時のフォールバック）
params.memory  = 16
params.threads = 4
params.target_region = 'v4'    // デフォルト（単体実行用）
params.output_layout = 'paired' 
params.read_length   = 250     

// 2. 固有の上限値定義
def BARRNAP_MAX_MEMORY  = 64
def BARRNAP_MAX_THREADS = 16

params.extract_16s_memory  = Math.min(params.memory as Integer, BARRNAP_MAX_MEMORY)
params.extract_16s_threads = Math.min(params.threads as Integer, BARRNAP_MAX_THREADS)

include { createNullParamsChannel; getParam; clusterOptions; processProfile; createPairsChannel; apptainerContainerOptions } \
    from "${params.petagenomeDir}/nf/common/utils"

process extract_16s_fastq {
    tag "${pair_id}_${region}"

    container = "${params.petagenomeDir}/modules/barrnap/barrnap.sif"
    containerOptions = { apptainerContainerOptions("${params.apptainerRunOptions}") }
    publishDir "${params.output}/${task.process}", mode: 'symlink', enabled: params.publish_output

    def gb = "${params.extract_16s_memory}"
    def threads = "${params.extract_16s_threads}"
    memory params.executor=="sge" ? null : "${gb} GB"
    cpus params.executor=="sge" ? null : threads
    clusterOptions "${clusterOptions(params.executor, gb, threads, label)}"

    input:
    // ★ pair_id, contigs に加えて region を入道させる
    tuple val(pair_id), path(contigs), val(region)

    output:
    // ★ 出力も領域ごとに独立したディレクトリ/ファイル名にする
    tuple val(pair_id), val(region), path("${pair_id}_${region}/${pair_id}_${region}*.fastq.gz")

    script:
    def layout = params.output_layout
    def read_len = params.read_length

    """
    echo "${processProfile(task)}" | tee prof.txt
    mkdir -p ${pair_id}_${region}

    # 1. barrnapで16S/23S等のrRNA領域を予測
    barrnap \
        --kingdom bac \
        --threads ${threads} \
        ${contigs} \
        > ${pair_id}_${region}/rRNA.gff

    # 2. Pythonスクリプトでインシリコ抽出 ＆ FASTQ変換
    python3 - << 'EOF'
import re
from Bio import SeqIO

region = "${region}"
layout = "${layout}"
read_len = ${read_len}
pair_id = "${pair_id}"
gff_file = "${pair_id}_${region}/rRNA.gff"
contig_file = "${contigs}"

PRIMERS = {
    'v1':  {'f': 'AGAGTTTGATCMTGGCTCAG', 'r': 'CTGCTGCCTCCCGTAGG'},        
    'v12': {'f': 'AGAGTTTGATCMTGGCTCAG', 'r': 'CGYCAATTCMTTTRWTTT'},   
    'v13': {'f': 'AGAGTTTGATCMTGGCTCAG', 'r': 'GWATTACCGCGGCKGCTG'},   
    'v34': {'f': 'CCTACGGGAGGCAGCAG',    'r': 'GGACTACNVGGGTWTCTAAT'}, 
    'v4':  {'f': 'GTGYCAGCMGCCGCGGTAA',    'r': 'GGACTACNVGGGTWTCTAAT'}  
}

if region not in PRIMERS:
    raise ValueError(f"Unknown target region: {region}. Choose from {list(PRIMERS.keys())}")

f_seq = PRIMERS[region]['f']
r_seq = PRIMERS[region]['r']

def iupac_to_regex(seq):
    mapping = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    return ''.join([mapping.get(b.upper(), 'N') for b in seq])

def rev_comp(seq):
    trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(trans)[::-1]

f_pat = re.compile(iupac_to_regex(f_seq), re.IGNORECASE)
r_pat = re.compile(iupac_to_regex(rev_comp(r_seq)), re.IGNORECASE)

contigs = SeqIO.to_dict(SeqIO.parse(contig_file, "fasta"))
amplicons = []

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or "16S" not in parts[8]: continue
        
        seqname, start, end, strand = parts[0], int(parts[3])-1, int(parts[4]), parts[6]
        if seqname not in contigs: continue
        
        seq = contigs[seqname].seq
        sub_seq = str(seq[start:end])
        if strand == "-":
            sub_seq = str(SeqIO.Seq(sub_seq).reverse_complement())
            
        for fm in f_pat.finditer(sub_seq):
            p1 = fm.start()
            search_win = sub_seq[p1:]
            for rm in r_pat.finditer(search_win):
                p2 = p1 + rm.end()
                amp = sub_seq[p1:p2]
                if len(amp) > 50: 
                    amplicons.append(amp)
                break

print(f"Extracted {len(amplicons)} in-silico amplicons for region {region}.")

import gzip

if layout == 'paired':
    r1_path = f"${pair_id}_${region}/${pair_id}_${region}_R1.fastq.gz"
    r2_path = f"${pair_id}_${region}/${pair_id}_${region}_R2.fastq.gz"
    
    with gzip.open(r1_path, 'wt') as f1, gzip.open(r2_path, 'wt') as f2:
        for i, amp in enumerate(amplicons):
            r1_seq = amp[:read_len]
            r2_seq = rev_comp(amp[-read_len:]) if len(amp) >= read_len else rev_comp(amp)
            qual = "I" * len(r1_seq)
            
            f1.write(f"@M03100:1:000000000-AMP:{i} 1:N:0:1\\n{r1_seq}\\n+\\n{qual}\\n")
            f2.write(f"@M03100:1:000000000-AMP:{i} 2:N:0:1\\n{r2_seq}\\n+\\n{qual}\\n")
else:
    r_path = f"${pair_id}_${region}/${pair_id}_${region}.fastq.gz"
    with gzip.open(r_path, 'wt') as f:
        for i, amp in enumerate(amplicons):
            qual = "I" * len(amp)
            f.write(f"@M03100:1:000000000-AMP:{i}\\n{amp}\\n+\\n{qual}\\n")

EOF
    """
}

workflow EXTRACT_16S_SUB {
    take:
    p
    contigs_with_region // [pair_id, contigs, region] のチャネルを想定

    main:
    out = extract_16s_fastq(contigs_with_region)

    emit:
    out = out
}

workflow EXTRACT_16S_ALL {
    p = createNullParamsChannel()
    contigs = createPairsChannel(params.extract_16s_contigs)
    regions_ch = Channel.from(params.target_regions.tokenize(',')).map { it.trim() }
    
    // 単体実行の際もコンティグと領域を組み合わせられるようにする
    combined_ch = contigs.combine(regions_ch)

    out_ch = EXTRACT_16S_SUB(p, combined_ch)
    out_ch.out.view { i -> "$i" }
}

workflow {
    EXTRACT_16S_ALL()
}