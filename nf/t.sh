#!/bin/bash

export TMPDIR=/dev/shm/tmp

threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=128
outdir=/dev/shm/petagenome_pipeline/out

data_dir="../modules/test"

ext_abs_dir="$(pwd)/../external"
virsorter_db="${ext_abs_dir}/virsorter-data"
virsorter_mga="${ext_abs_dir}/mga_linux_ia64"
virsorter2_db="${ext_abs_dir}/virsorter2-data"
metaphlan2_db="${ext_abs_dir}/metaphlan2_db"
metaphlan4_db="${ext_abs_dir}/metaphlan4_db"

args="\
    --petagenomeDir=$(pwd)/.. \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

test=${1:-"main"}

case ${test} in
    "main")
        nextflow run main.nf ${args} \
                 --test_main_reads "${data_dir}/ecoli_1K_{1,2}.fq.gz"
        ;;
    "bbmap")
        nextflow run bbmap.nf ${args} \
                 --test_bbmap_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_bbmap_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "blast")
        nextflow run blast.nf ${args} \
                 --test_blast_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_blast_qry "${data_dir}/q1.fasta"
        ;;
    "bowtie")
        nextflow run bowtie.nf ${args} \
                 --test_bowtie_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_bowtie_qry "${data_dir}/q1.fasta"
        ;;
    "bowtie2")
        nextflow run bowtie2.nf ${args} \
                 --test_bowtie2_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_bowtie2_qry "${data_dir}/q1.fasta"
        ;;
    "bwa")
        nextflow run bwa.nf ${args} \
                 --test_bwa_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_bwa_qry "${data_dir}/q1.fasta"
        ;;
    "bwa_mem2")
        nextflow run bwa_mem2.nf ${args} \
                 --test_bwa_mem2_ref "${data_dir}/ecoli_1K_1.fa.gz" \
                 --test_bwa_mem2_qry "${data_dir}/q1.fasta"
        ;;
    "cdhit")
        nextflow run cdhit.nf ${args} \
                 --test_cdhit_read "${data_dir}/ecoli_1K_1.fa.gz"
        ;;
    "cutadapt")
        nextflow run cutadapt.nf ${args} \
                 --test_cutadapt_reads "${data_dir}/ecoli_1K_{1,2}.fq.gz"
        ;;
    "diamond")
        nextflow run diamond.nf ${args} \
                 --test_diamond_ref "${data_dir}/1.faa.gz" \
                 --test_diamond_qry "${data_dir}/2.faa"
        ;;
    "falco")
        nextflow run falco.nf ${args} \
                 --test_falco_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "fastp")
        nextflow run fastp.nf ${args} \
                 --test_fastp_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "fastqc")
        nextflow run fastqc.nf ${args} \
                 --test_fastqc_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "megahit")
        nextflow run megahit.nf ${args} \
                 --test_megahit_reads "${data_dir}/ecoli_1K_{1,2}.fq.gz"
        ;;
    "metaphlan")
        nextflow run metaphlan.nf ${args} \
                 --test_metaphlan_read "${data_dir}/s_6_1.fastq.gz" \
                 --metaphlan_db ${metaphlan4_db}
        ;;
    "metaphlan2")
        nextflow run metaphlan2.nf ${args} \
                 --test_metaphlan2_read "${data_dir}/s_6_1.fastq.gz" \
                 --metaphlan2_db ${metaphlan2_db}
        ;;
    "minimap2")
        nextflow run minimap2.nf ${args} \
                 --test_minimap2_ref "${data_dir}/8seq.fa" \
                 --test_minimap2_qry "${data_dir}/1seq.fa"
        ;;
    "prinseq")
        nextflow run prinseq.nf ${args} \
                 --test_prinseq_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "prodigal")
        nextflow run prodigal.nf ${args} \
                 --test_prodigal_read "${data_dir}/ecoli_1K_1.fa.gz"
        ;;
    "spades")
        nextflow run spades.nf ${args} \
                 --test_spades_reads "${data_dir}/ecoli_1K_{1,2}.fq.gz"
        nextflow run spades.nf \
                 --test_spades_reads "${data_dir}/s_6_{1,2}.fastq.gz"
        ;;
    "virsorter")
        nextflow run virsorter.nf ${args} \
                 --virsorter_db ${virsorter_db} \
                 --virsorter_mga ${virsorter_mga} \
                 --virsorter_db_type refseq \
                 --virsorter_aligner blast
                 --test_virsorter_read "${data_dir}/1seq.fa" \
        nextflow run virsorter.nf ${args} \
                 --virsorter_db ${virsorter_db} \
                 --virsorter_mga ${virsorter_mga} \
                 --virsorter_db_type virome \
                 --virsorter_aligner diamond
                 --test_virsorter_read "${data_dir}/1seq.fa" \
        ;;
    "virsorter2")
        nextflow run virsorter2.nf ${args} \
                 --virsorter2_db ${virsorter2_db} \
                 --test_virsorter2_read "${data_dir}/1seq.fa"
        ;;
    "*")
esac
