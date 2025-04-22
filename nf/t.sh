#!/bin/bash

#nextflow run main.nf --reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run fastp.nf --test_fastp_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run spades.nf --test_spades_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run cutadapt.nf --test_cutadapt_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --petagenomeDir=$(pwd)/..

#nextflow run prinseq.nf --test_prinseq_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..

#nextflow run bbmap.nf --test_bbmap_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bbmap_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..

#nextflow run cdhit.nf --test_cdhit_read "../modules/test/ecoli_1K_1.fa.gz" --petagenomeDir=$(pwd)/..

#nextflow run megahit.nf --test_megahit_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run blast.nf --test_blast_ref "../modules/test/ecoli_1K_1.fa.gz" --test_blast_qry "../modules/test/q1.fasta" --petagenomeDir=$(pwd)/..

nextflow run diamond.nf --test_diamond_ref "../modules/test/1.faa" --test_diamond_qry "../modules/test/2.faa" --petagenomeDir=$(pwd)/..

################################# NG
#nextflow run spades.nf --test_spades_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

