#!/bin/bash

export TMPDIR=/dev/shm/tmp

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0

nextflow run main.nf --output "out" --reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..



#nextflow run bbmap.nf --test_bbmap_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bbmap_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..

#nextflow run blast.nf --test_blast_ref "../modules/test/ecoli_1K_1.fa.gz" --test_blast_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run bowtie.nf --test_bowtie_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bowtie_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --random_seed ${random_seed} --petagenomeDir=$(pwd)/..

#nextflow run bowtie2.nf --test_bowtie2_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bowtie2_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --random_seed ${random_seed} --petagenomeDir=$(pwd)/..

#nextflow run bwa.nf --test_bwa_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bwa_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run bwa_mem2.nf --test_bwa_mem2_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bwa_mem2_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run cdhit.nf --test_cdhit_read "../modules/test/ecoli_1K_1.fa.gz" --petagenomeDir=$(pwd)/..

#nextflow run cutadapt.nf --test_cutadapt_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --petagenomeDir=$(pwd)/..

#nextflow run diamond.nf --test_diamond_ref "../modules/test/1.faa.gz" --test_diamond_qry "../modules/test/2.faa" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run falco.nf --test_falco_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run fastp.nf --test_fastp_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run fastqc.nf --test_fastqc_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run megahit.nf --test_megahit_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run metaphlan.nf --test_metaphlan_read "../modules/test/s_6_1.fastq.gz" --metaphlan_db "$(pwd)/../external/metaphlan4_db" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run metaphlan2.nf --test_metaphlan2_read "../modules/test/s_6_1.fastq.gz" --metaphlan2_db "$(pwd)/../external/metaphlan2_db" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run minimap2.nf --test_minimap2_ref "../modules/test/8seq.fa" --test_minimap2_qry "../modules/test/1seq.fa" --petagenomeDir=$(pwd)/..

#nextflow run prinseq.nf --test_prinseq_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..

#nextflow run prodigal.nf --test_prodigal_read "../modules/test/ecoli_1K_1.fa.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run spades.nf --test_spades_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run spades.nf --test_spades_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#nextflow run virsorter.nf --virsorter_db "$(pwd)/../external/virsorter-data" --virsorter_mga "$(pwd)/../external/mga_linux_ia64" --test_virsorter_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --virsorter_db_type refseq --virsorter_aligner blast --petagenomeDir=$(pwd)/..

#nextflow run virsorter.nf --virsorter_db "$(pwd)/../external/virsorter-data" --virsorter_mga "$(pwd)/../external/mga_linux_ia64" --test_virsorter_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --virsorter_db_type virome --virsorter_aligner diamond --petagenomeDir=$(pwd)/..

#nextflow run virsorter2.nf --virsorter2_db "$(pwd)/../external/virsorter2-data" --test_virsorter2_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..

#------

#



################################# NG

