#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ncbi-blast.sif ncbi-blast -h

exit

# #{blast_dir}/makeblastdb -in #{out}.1000.fa -out #{out}.1000 -dbtype nucl -parse_seqids
# #{blast_dir}/bin/blastn -task #{task} -num_threads #{n_threads} -query #{prophage}.fa -db #{contig_proph} -perc_identity #{pi_thre} -evalue #{e_thre} -outfmt 6 -num_alignments #{n_al} -out #{blast_out}.tx
# 03_assembly.rb 04_pool_contigs.rb 05_circular_contigs.rb 06_prophage_detection.02_after_processing.rb











