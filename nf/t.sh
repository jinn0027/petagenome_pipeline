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

nfDir="${PETAGENOME_PIPELINE_DIR}/nf"
dataDir="${PETAGENOME_PIPELINE_DIR}/modules/test"
extDir="${PETAGENOME_PIPELINE_DIR}/external"

virsorterDb="${extDir}/virsorter-data"
virsorterMga="${extDir}/mga_linux_ia64"
virsorter2Db="${extDir}/virsorter2-data"
metaphlanDb="${extDir}/metaphlan_db"

longFnqGzPair="${dataDir}/ecoli_1K_{1,2}.fq.gz"
longFnaPair="${dataDir}/NC_*.fna"
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

args+=" --publish_output true"

#test=${1:-"error_correction"}
#test=${1:-"assembly"}
#test=${1:-"pool_contigs"}
test=${1:-"circular_contigs"}

test=${test%.*}

case ${test} in
    "main")
        nextflow run ${nfDir}/toys/main.nf ${args} \
                 --test_main_reads "${longFnqGzPair}"
        ;;
    "error_correction")
        nextflow run ${nfDir}/lv2/error_correction.nf ${args} \
                 --test_error_correction_reads "${longFnqGzPair}"
        ;;
    "assembly")
        nextflow run ${nfDir}/lv2/assembly.nf ${args} \
                 --test_assembly_reads "${longFnqGzPair}"
        ;;
    "pool_contigs")
        nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
                 -with-report report_${test}.html \
                 -with-trace trace_${test}.txt \
                 --test_pool_contigs_contigs "${longFnaPair}"
        ;;
    "circular_contigs")
        nextflow run ${nfDir}/lv2/circular_contigs.nf ${args} \
                 --test_circular_contigs_contig "${longFnaGz1}"
        ;;
    "bbmap")
        nextflow run ${nfDir}/lv1/bbmap.nf ${args} \
                 --test_bbmap_ref "${longFnaGz1}" \
                 --test_bbmap_reads "${shortFnqGzPair}"
        ;;
    "blast")
        nextflow run ${nfDir}/lv1/blast.nf ${args} \
                 --test_blast_ref "${longFnaGz1}" \
                 --test_blast_qry "${shortFna1}"
        ;;
    "bowtie")
        nextflow run ${nfDir}/lv1/bowtie.nf ${args} \
                 --test_bowtie_ref "${longFnaGz1}" \
                 --test_bowtie_qry "${shortFna1}"
        ;;
    "bowtie2")
        nextflow run ${nfDir}/lv1/bowtie2.nf ${args} \
                 --test_bowtie2_ref "${longFnaGz1}" \
                 --test_bowtie2_qry "${shortFna1}"
        ;;
    "bwa")
        nextflow run ${nfDir}/lv1/bwa.nf ${args} \
                 --test_bwa_ref "${longFnaGz1}" \
                 --test_bwa_qry "${shortFna1}"
        ;;
    "bwa_mem2")
        nextflow run ${nfDir}/lv1/bwa_mem2.nf ${args} \
                 --test_bwa_mem2_ref "${longFnaGz1}" \
                 --test_bwa_mem2_qry "${shortFna1}"
        ;;
    "cdhit")
        nextflow run ${nfDir}/lv1/cdhit.nf ${args} \
                 --test_cdhit_read "${shortFnaGz1}"
        ;;
    "cutadapt")
        nextflow run ${nfDir}/lv1/cutadapt.nf ${args} \
                 --test_cutadapt_reads "${longFnqGzPair}"
        ;;
    "diamond")
        nextflow run ${nfDir}/lv1/diamond.nf ${args} \
                 --test_diamond_ref "${shortFaaGz1}" \
                 --test_diamond_qry "${shortFaa2}"
        ;;
    "falco")
        nextflow run ${nfDir}/lv1/falco.nf ${args} \
                 --test_falco_reads "${shortFnqGzPair}"
        ;;
    "fastp")
        nextflow run ${nfDir}/lv1/fastp.nf ${args} \
                 --test_fastp_reads "${shortFnqGzPair}"
        ;;
    "fastqc")
        nextflow run ${nfDir}/lv1/fastqc.nf ${args} \
                 --test_fastqc_reads "${shortFnqGzPair}"
        ;;
    "megahit")
        nextflow run ${nfDir}/lv1/megahit.nf ${args} \
                 --test_megahit_reads "${longFnqGzPair}"
        ;;
    "metaphlan")
        nextflow run ${nfDir}/lv1/metaphlan.nf ${args} \
                 --test_metaphlan_read "${shortFnqGz1}" \
                 --metaphlan_db "${metaphlanDb}"
        ;;
    "minimap2")
        nextflow run ${nfDir}/lv1/minimap2.nf ${args} \
                 --test_minimap2_ref "${longFnaX8}" \
                 --test_minimap2_qry "${longFnaX1}"
        ;;
    "mmseqs2")
        nextflow run ${nfDir}/lv1/mmseqs2.nf ${args} \
                 --test_mmseqs2_ref "${longFnaX8}" \
                 --test_mmseqs2_qry "${longFnaX1}"
        ;;
    "prinseq")
        nextflow run ${nfDir}/lv1/prinseq.nf ${args} \
                 --test_prinseq_reads "${shortFnqGzPair}"
        ;;
    "prodigal")
        nextflow run ${nfDir}/lv1/prodigal.nf ${args} \
                 --test_prodigal_read "${longFnaGz1}"
        ;;
    "spades")
        nextflow run ${nfDir}/lv1/spades.nf ${args} \
                 --test_spades_reads "${longFnqGzPair}"
        nextflow run ${nfDir}/lv1/spades.nf ${args} \
                 --test_spades_reads "${shortFnqGzPair}"
        ;;
    "virsorter")
        nextflow run ${nfDir}/lv1/virsorter.nf ${args} \
                 --virsorter_db "${virsorterDb}" \
                 --virsorter_mga "${virsorterMga}" \
                 --virsorter_db_type refseq \
                 --virsorter_aligner blast \
                 --test_virsorter_read "${longFnaX1}" \
        nextflow run ${nfDir}/lv1/virsorter.nf ${args} \
                 --virsorter_db "${virsorterDb}" \
                 --virsorter_mga "${virsorterMga}" \
                 --virsorter_db_type virome \
                 --virsorter_aligner diamond \
                 --test_virsorter_read "${longFnaX1}" \
        ;;
    "virsorter2")
        nextflow run ${nfDir}/lv1/virsorter2.nf ${args} \
                 --virsorter2_db "${virsorter2Db}" \
                 --test_virsorter2_read "${longFnaX1}"
        ;;
    "helloruby")
        nextflow run ${nfDir}/toys/helloruby.nf ${args} \
                 --test_helloruby_reads "${shortFnqGzPair}"
        ;;
    "*")
esac

#rm -rf /dev/shm/${USER}
