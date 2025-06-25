#!/bin/bash

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

threads=16
#threads=$(nproc)
cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
random_seed=0
memory=32
#memory=128
#outdir=/dev/shm/${USER}/petagenome_pipeline/out
outdir=out

#dataDir="${MY_DIR}/../modules/test"
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"
#extDir="${MY_DIR}/../external"
extDir="${PETAGENOME_PIPELINE_DIR}/external"

virsorterDb="${extDir}/virsorter-data"
virsorterMga="${extDir}/mga_linux_ia64"
virsorter2Db="${extDir}/virsorter2-data"
metaphlan2Db="${extDir}/metaphlan2_db"
metaphlan4Db="${extDir}/metaphlan4_db"

longFnqGzPair="${dataDir}/ecoli_1K_{1,2}.fq.gz"
longFnaGz1="${dataDir}/ecoli_1K_1.fa.gz"
longFnaGz2="${dataDir}/ecoli_1K_2.fa.gz"
longFnaX1="${dataDir}/1seq.fa"
longFnaX8="${dataDir}/8seq.fa"
shortFnqGzPair="${dataDir}/s_6_{1,2}.fastq.gz"
shortFnaGz1="${dataDir}/s_6_1.fasta.gz"
shortFnqGz1="${dataDir}/s_6_1.fastq.gz"
shortFaaGz1="${dataDir}/1.faa.gz"
shortFaa2="${dataDir}/2.faa"
shortFna1="${dataDir}/q.fa"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

test=${1:-"error_correction"}
#test=${1:-"assembly"}

case ${test} in
    "main")
        nextflow run main.nf ${args} \
                 --test_main_reads ${longFnqGzPair}
        ;;
    "error_correction")
        nextflow run error_correction.nf ${args} \
                 --assembly_type virome \
                 --test_error_correction_reads ${longFnqGzPair}
        ;;
    "assembly")
        nextflow run assembly.nf ${args} \
                 --assembly_type virome \
                 --test_assembly_reads ${longFnqGzPair}
        ;;
    "bbmap")
        nextflow run bbmap.nf ${args} \
                 --test_bbmap_ref ${longFnaGz1} \
                 --test_bbmap_reads ${shortFnqGzPair}
        ;;
    "blast")
        nextflow run blast.nf ${args} \
                 --test_blast_ref ${longFnaGz1} \
                 --test_blast_qry ${shortFna1}
        ;;
    "bowtie")
        nextflow run bowtie.nf ${args} \
                 --test_bowtie_ref ${longFnaGz1} \
                 --test_bowtie_qry ${shortFna1}
        ;;
    "bowtie2")
        nextflow run bowtie2.nf ${args} \
                 --test_bowtie2_ref ${longFnaGz1} \
                 --test_bowtie2_qry ${shortFna1}
        ;;
    "bwa")
        nextflow run bwa.nf ${args} \
                 --test_bwa_ref ${longFnaGz1} \
                 --test_bwa_qry ${shortFna1}
        ;;
    "bwa_mem2")
        nextflow run bwa_mem2.nf ${args} \
                 --test_bwa_mem2_ref ${longFnaGz1} \
                 --test_bwa_mem2_qry ${shortFna1}
        ;;
    "cdhit")
        nextflow run cdhit.nf ${args} \
                 --test_cdhit_read ${shortFnaGz1}
        ;;
    "cutadapt")
        nextflow run cutadapt.nf ${args} \
                 --test_cutadapt_reads ${longFnqGzPair}
        ;;
    "diamond")
        nextflow run diamond.nf ${args} \
                 --test_diamond_ref ${shortFaaGz1} \
                 --test_diamond_qry ${shortFaa2}
        ;;
    "falco")
        nextflow run falco.nf ${args} \
                 --test_falco_reads ${shortFnqGzPair}
        ;;
    "fastp")
        nextflow run fastp.nf ${args} \
                 --test_fastp_reads ${shortFnqGzPair}
        ;;
    "fastqc")
        nextflow run fastqc.nf ${args} \
                 --test_fastqc_reads ${shortFnqGzPair}
        ;;
    "megahit")
        nextflow run megahit.nf ${args} \
                 --test_megahit_reads ${longFnqGzPair}
        ;;
    "metaphlan")
        nextflow run metaphlan.nf ${args} \
                 --test_metaphlan_read ${shortFnqGz1} \
                 --metaphlan_db ${metaphlan4Db}
        ;;
    "metaphlan2")
        nextflow run metaphlan2.nf ${args} \
                 --test_metaphlan2_read ${shortFnqGz1} \
                 --metaphlan2_db ${metaphlan2Db}
        ;;
    "minimap2")
        nextflow run minimap2.nf ${args} \
                 --test_minimap2_ref ${longFnaX8} \
                 --test_minimap2_qry ${longFnaX1}
        ;;
    "mmseqs2")
        nextflow run mmseqs2.nf ${args} \
                 --test_mmseqs2_ref ${longFnaX8} \
                 --test_mmseqs2_qry ${longFnaX1}
        ;;
    "prinseq")
        nextflow run prinseq.nf ${args} \
                 --test_prinseq_reads ${shortFnqGzPair}
        ;;
    "prodigal")
        nextflow run prodigal.nf ${args} \
                 --test_prodigal_read ${longFnaGz1}
        ;;
    "spades")
        nextflow run spades.nf ${args} \
                 --test_spades_reads ${longFnqGzPair}
        nextflow run spades.nf \
                 --test_spades_reads ${shortFnqGzPair}
        ;;
    "virsorter")
        nextflow run virsorter.nf ${args} \
                 --virsorter_db ${virsorterDb} \
                 --virsorter_mga ${virsorterMga} \
                 --virsorter_db_type refseq \
                 --virsorter_aligner blast \
                 --test_virsorter_read ${longFnaX1} \
        nextflow run virsorter.nf ${args} \
                 --virsorter_db ${virsorterDb} \
                 --virsorter_mga ${virsorterMga} \
                 --virsorter_db_type virome \
                 --virsorter_aligner diamond \
                 --test_virsorter_read ${longFnaX1} \
        ;;
    "virsorter2")
        nextflow run virsorter2.nf ${args} \
                 --virsorter2_db ${virsorter2Db} \
                 --test_virsorter2_read ${longFnaX1}
        ;;
    "toys/helloruby")
        nextflow run toys/helloruby.nf ${args} \
                 --test_helloruby_reads ${shortFnqGzPair}
        ;;
    "*")
esac

rm -rf /dev/shm/${USER}
