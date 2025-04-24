#!/bin/bash

export TMPDIR=/dev/shm/tmp

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
outdir=/dev/shm/petagenome_pipeline/out

test=${1:-"main"}

case ${test} in
    "main")
        nextflow run main.nf --output ${outdir} --reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "bbmap")
        nextflow run bbmap.nf --output ${outdir} --test_bbmap_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bbmap_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..
        ;;
    "blast")
        nextflow run blast.nf --output ${outdir} --test_blast_ref "../modules/test/ecoli_1K_1.fa.gz" --test_blast_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "bowtie")
        nextflow run bowtie.nf --output ${outdir} --test_bowtie_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bowtie_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --random_seed ${random_seed} --petagenomeDir=$(pwd)/..
        ;;
    "bowtie2")
        nextflow run bowtie2.nf --output ${outdir} --test_bowtie2_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bowtie2_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --random_seed ${random_seed} --petagenomeDir=$(pwd)/..
        ;;
    "bwa")
        nextflow run bwa.nf --output ${outdir} --test_bwa_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bwa_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "bwa_mem2")
        nextflow run bwa_mem2.nf --output ${outdir} --test_bwa_mem2_ref "../modules/test/ecoli_1K_1.fa.gz" --test_bwa_mem2_qry "../modules/test/q1.fasta" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "cdhit")
        nextflow run cdhit.nf --output ${outdir} --test_cdhit_read "../modules/test/ecoli_1K_1.fa.gz" --petagenomeDir=$(pwd)/..
        ;;
    "cutadapt")
        nextflow run cutadapt.nf --output ${outdir} --test_cutadapt_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --petagenomeDir=$(pwd)/..
        ;;
    "diamond")
        nextflow run diamond.nf --output ${outdir} --test_diamond_ref "../modules/test/1.faa.gz" --test_diamond_qry "../modules/test/2.faa" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "falco")
        nextflow run falco.nf --output ${outdir} --test_falco_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "fastp")
        nextflow run fastp.nf --output ${outdir} --test_fastp_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "fastqc")
        nextflow run fastqc.nf --output ${outdir} --test_fastqc_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "megahit")
        nextflow run megahit.nf --output ${outdir} --test_megahit_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "metaphlan")
        nextflow run metaphlan.nf --output ${outdir} --test_metaphlan_read "../modules/test/s_6_1.fastq.gz" --metaphlan_db "$(pwd)/../external/metaphlan4_db" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "metaphlan2")
        nextflow run metaphlan2.nf --output ${outdir} --test_metaphlan2_read "../modules/test/s_6_1.fastq.gz" --metaphlan2_db "$(pwd)/../external/metaphlan2_db" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "minimap2")
        nextflow run minimap2.nf --output ${outdir} --test_minimap2_ref "../modules/test/8seq.fa" --test_minimap2_qry "../modules/test/1seq.fa" --petagenomeDir=$(pwd)/..
        ;;
    "prinseq")
        nextflow run prinseq.nf --output ${outdir} --test_prinseq_reads "../modules/test/s_6_{1,2}.fastq.gz" --petagenomeDir=$(pwd)/..
        ;;
    "prodigal")
        nextflow run prodigal.nf --output ${outdir} --test_prodigal_read "../modules/test/ecoli_1K_1.fa.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "spades")
        nextflow run spades.nf --output ${outdir} --test_spades_reads "../modules/test/ecoli_1K_{1,2}.fq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        nextflow run spades.nf --output ${outdir} --test_spades_reads "../modules/test/s_6_{1,2}.fastq.gz" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "virsorter")
        nextflow run virsorter.nf --output ${outdir} --virsorter_db "$(pwd)/../external/virsorter-data" --virsorter_mga "$(pwd)/../external/mga_linux_ia64" --test_virsorter_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --virsorter_db_type refseq --virsorter_aligner blast --petagenomeDir=$(pwd)/..
        nextflow run virsorter.nf --output ${outdir} --virsorter_db "$(pwd)/../external/virsorter-data" --virsorter_mga "$(pwd)/../external/mga_linux_ia64" --test_virsorter_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --virsorter_db_type virome --virsorter_aligner diamond --petagenomeDir=$(pwd)/..
        ;;
    "virsorter2")
        nextflow run virsorter2.nf --output ${outdir} --virsorter2_db "$(pwd)/../external/virsorter2-data" --test_virsorter2_read "../modules/test/1seq.fa" --threads ${threads} --cpus ${cpus} --petagenomeDir=$(pwd)/..
        ;;
    "*")
esac
