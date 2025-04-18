#!/bin/bash

#nextflow run main.nf --reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run fastp.nf --test_fastp_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

#nextflow run spades.nf --test_spades_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

nextflow run cutadapt.nf --test_cutadapt_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --petagenomeDir=$(pwd)/..
