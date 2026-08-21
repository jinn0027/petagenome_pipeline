#!/bin/bash

if [[ ! -v PETAGENOME_PIPELINE_DIR ]] ; then
    echo "PETAGENOME_PIPELINE_DIR not defined"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
if [ ! -d ${PETAGENOME_PIPELINE_DIR} ] ; then
    echo "${PETAGENOME_PIPELINE_DIR} does not exist"
    echo "Please source <petagenome_dir>/etc/host_setup.sh"
    exit 1
fi
echo "PETAGENOME_PIPELINE_DIR : ${PETAGENOME_PIPELINE_DIR}"

#export TMPDIR=/dev/shm/${USER}/tmp
export TMPDIR=$(pwd)/tmp

mkdir -p ${TMPDIR}

MY_FILE="${BASH_SOURCE[0]}"
MY_DIR="$(cd "$(dirname "${MY_FILE}")" && pwd)"

#date=$(date +"%Y%m%d%H%M%S")
date="aaa"

#threads=16
threads=64
#threads=$(nproc)
#cpus=$(grep physical.id /proc/cpuinfo | sort -u | wc -l)
cpus=16
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
shortFnaX2Gz1="${dataDir}/s_6_1_x2.fasta.gz"
shortFaaGz1="${dataDir}/1.faa.gz"
shortFaa2="${dataDir}/2.faa"
shortFna1="${dataDir}/q.fa"
aFna="${dataDir}/a.fna"
aFaa="${dataDir}/a.faa"
bFna="${dataDir}/b.fna"
bFaa="${dataDir}/b.faa"
cFna="${dataDir}/c.fna"
cFaa="${dataDir}/c.faa"
dFna="${dataDir}/d.fna"
dFaa="${dataDir}/d.faa"

inPairs="/home/data202502/metagenome/*_XXXXXXXX_XXXXXXXX_L001_R{1,2}_001.fastq.gz"

args="\
    --petagenomeDir=${PETAGENOME_PIPELINE_DIR} \
    --output ${outdir} \
    --memory ${memory} \
    --threads ${threads} \
    --cpus ${cpus} \
    --random_seed ${random_seed} \
    "

args+=" --publish_output true"
args+=" -profile slurm"

args_dbg="\
    -with-trace trace.${date}.txt \
    -with-report report.${date}.html \
    -with-timeline timeline.${date}.html \
    -with-dag dag.${date}.png \
    "

args+=" ${args_dbg}"

test=${1:-"remove_host"}
#test=${1:-"aa"}
#test=${1:-"error_correction"}
#test=${1:-"assembly"}
#test=${1:-"pool_contigs"}
#test=${1:-"circular_contigs"}

#nextflow clean -f

test=${test%.*}

case ${test} in
    "fastp")
        nextflow run ${nfDir}/lv1/fastp.nf -entry FASTP_ALL ${args} \
                 --fastp_reads "${shortFnqGzPair}"
        ;;
    "remove_host")
        db=${extDir}/GRCh38/bwa_db
        if [ ! -d ${db} ] ; then
            ref=${extDir}/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
            if [ ! -f ${ref} ] ; then
                echo "Error : ${ref} not found."
                echo "Please download as:"
                echo "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
                exit 1
            fi
            nextflow run ${nfDir}/lv2/remove_host.nf -entry BUILD_REF_DB_ONLY ${args} \
                     --remove_host_aligner "bwa_mem2" \
                     --remove_host_ref_fasta "${ref}"
            mv ${outdir}/BUILD_REF_DB_ONLY:BUILD_REF_DB_SUB:bwa_mem2_makerefdb/00_GCF ${db}
        fi
        nextflow run ${nfDir}/lv2/remove_host.nf -entry REMOVE_HOST_WITH_DB ${args} \
                 --remove_host_aligner "bwa_mem2" \
                 --remove_host_is_prebuilt_db "true" \
                 --remove_host_prebuilt_db ${db} \
                 --remove_host_reads "${shortFnqGzPair}"
        ;;
    "assembly")
        nextflow run ${nfDir}/lv2/assembly.nf -entry ASSEMBLY_ALL ${args} \
                 --assembly_assembler "megahit" \
                 --assembly_l_thre 10 \
                 --assembly_reads "${longFnqGzPair}"
        ;;
    "prodigal")
        nextflow run ${nfDir}/lv1/prodigal.nf -entry PRODIGAL_ALL ${args} \
                 --prodigal_procedure "meta" \
                 --prodigal_read "${longFnaGz1}"
        ;;
    "annotate_catalog_prot")
        db=${extDir}/uniprot_refs/mmseqs2
        if [ ! -d ${db} ] ; then
            mkdir -p $(dirname ${db})
            ref=${extDir}/uniprot_refs/uniprot_sprot.fasta.gz 
            nextflow run ${nfDir}/lv2/annotate_catalog_prot.nf -entry BUILD_REF_DB_ONLY ${args} \
                     --annotate_catalog_prot_aligner mmseqs2 \
                     --annotate_catalog_prot_ref_fasta ${ref} \
                     --mmseqs2_ref_type 1
            mv ${outdir}/BUILD_REF_DB_ONLY:BUILD_REF_DB_SUB:mmseqs2_makerefdb/00_uniprot ${db}
        fi
        nextflow run ${nfDir}/lv2/annotate_catalog_prot.nf -entry ANNOTATE_WITH_DB ${args} \
                 --annotate_catalog_prot_aligner mmseqs2 \
                 --annotate_catalog_prot_is_prebuilt_db "true" \
                 --annotate_catalog_prot_prebuilt_db ${db} \
                 --annotate_catalog_prot_orfs ${shortFnaGz1} \
                 --annotate_catalog_prot_taxid_map ${extDir}/uniprot_refs/uniprot_to_taxid.tsv \
                 --annotate_catalog_prot_ko_map ${extDir}/uniprot_refs/uniprot_to_ko.tsv \
                 --prodigal_procedure "meta" \
                 --mmseqs2_ref_type 1 \
                 --mmseqs2_qry_type 1
        ;;
    "nr_catalog")
        nextflow run ${nfDir}/lv2/nr_catalog.nf -entry NR_CATALOG_ALL ${args} \
                 --nr_catalog_faa "${aFaa};${bFaa};${cFaa};${dFaa}" \
                 --nr_catalog_fna "${aFna};${bFna};${cFna};${dFna}"
        ;;
    "aa")
        nextflow run ${nfDir}/lv3/aa.nf -resume ${args} \
                 --remove_host_aligner "bwa_mem2" \
                 --remove_host_is_prebuilt_db "true" \
                 --remove_host_ref_fasta_or_db "${extDir}/GRCh38/bwa_db" \
                 --assembly_assembler "megahit" \
                 --annotate_catalog_prot_aligner mmseqs2 \
                 --annotate_catalog_prot_is_prebuilt_db "true" \
                 --annotate_catalog_prot_ref_or_db ${extDir}/uniprot_refs/mmseqs2 \
                 --taxid_map_path ${extDir}/uniprot_refs/uniprot_to_taxid.tsv \
                 --ko_map_path ${extDir}/uniprot_refs/uniprot_to_ko.tsv \
                 --taxid_name_map_path ${extDir}/uniprot_refs/taxid_to_name.tsv \
                 --ko_name_map_path ${extDir}/uniprot_refs/ko_to_name.tsv \
                 --bacteriome_pipeline_reads "${inPairs}"

#                 --assembly_l_thre 10 \
#                 --mmseqs2_ref_type 1 \
#                 --mmseqs2_qry_type 1 \
#                 --host_is_prebuilt_db "false" \
#                 --bacteriome_pipeline_reads "${longFnqGzPair}"
#                 --host_ref_fasta_or_db "${dataDir}/ecoli_1K_1.fa.gz" \
        ;;
    "main")
        nextflow run ${nfDir}/toys/main.nf ${args} \
                 --main_reads "${longFnqGzPair}"
        ;;
    "error_correction")
        nextflow run ${nfDir}/lv2/error_correction.nf ${args} \
                 --error_correction_reads "${longFnqGzPair}"
        ;;
    "pool_contigs")
        nextflow run ${nfDir}/lv2/pool_contigs.nf ${args} \
                 --pool_contigs_l_thre 10 \
                 --pool_contigs_contigs "${longFnaPair}"
        ;;
    "circular_contigs")
        nextflow run ${nfDir}/lv2/circular_contigs.nf ${args} \
                 --circular_contigs_contig "${longFnaGz1}"
        ;;
    "bbmap")
        nextflow run ${nfDir}/lv1/bbmap.nf ${args} \
                 --bbmap_ref "${longFnaGz1}" \
                 --bbmap_reads "${shortFnqGzPair}"
        ;;
    "blast")
        #nextflow run ${nfDir}/lv1/blast.nf -entry BUILD_REF_DB_ONLY ${args} \
            nextflow run ${nfDir}/lv1/blast.nf ${args} \
                 --blast_ref "${longFnaGz1}" \
                 --blast_qry "${shortFna1}"
        ;;
    "bowtie")
        nextflow run ${nfDir}/lv1/bowtie.nf ${args} \
                 --bowtie_ref "${longFnaGz1}" \
                 --bowtie_qry "${shortFna1}"
        ;;
    "bowtie2")
        nextflow run ${nfDir}/lv1/bowtie2.nf ${args} \
                 --bowtie2_ref "${longFnaGz1}" \
                 --bowtie2_qry "${shortFna1}"
        ;;
    "bwa")
        nextflow run ${nfDir}/lv1/bwa.nf ${args} \
                 --bwa_ref "${longFnaGz1}" \
                 --bwa_qry "${shortFna1}"
        ;;
    "bwa_mem2")
        nextflow run ${nfDir}/lv1/bwa_mem2.nf ${args} \
                 --bwa_mem2_ref "${longFnaGz1}" \
                 --bwa_mem2_qry "${shortFna1}"
        ;;
    "cdhit")
        nextflow run ${nfDir}/lv1/cdhit.nf ${args} \
                 --cdhit_read "${shortFnaX2Gz1}"
        ;;
    "cutadapt")
        nextflow run ${nfDir}/lv1/cutadapt.nf ${args} \
                 --cutadapt_reads "${longFnqGzPair}"
        ;;
    "diamond")
        nextflow run ${nfDir}/lv1/diamond.nf ${args} \
                 --diamond_ref "${shortFaaGz1}" \
                 --diamond_qry "${shortFaa2}"
        ;;
    "falco")
        nextflow run ${nfDir}/lv1/falco.nf ${args} \
                 --falco_reads "${shortFnqGzPair}"
        ;;
    "fastqc")
        nextflow run ${nfDir}/lv1/fastqc.nf ${args} \
                 --fastqc_reads "${shortFnqGzPair}"
        ;;
    "megahit")
        nextflow run ${nfDir}/lv1/megahit.nf ${args} \
                 --megahit_reads "${longFnqGzPair}"
        ;;
    "metaphlan")
        nextflow run ${nfDir}/lv1/metaphlan.nf ${args} \
                 --metaphlan_read "${shortFnqGz1}" \
                 --metaphlan_db "${metaphlanDb}"
        ;;
    "minimap2")
        nextflow run ${nfDir}/lv1/minimap2.nf ${args} \
                 --minimap2_ref "${longFnaX8}" \
                 --minimap2_qry "${longFnaX1}"
        ;;
    "mmseqs2")
        nextflow run ${nfDir}/lv1/mmseqs2.nf ${args} \
                 --mmseqs2_cluster_mode cluster \
                 --mmseqs2_module cluster \
                 --mmseqs2_ref "${shortFnaX2Gz1}"
        nextflow run ${nfDir}/lv1/mmseqs2.nf ${args} \
                 --mmseqs2_search_type 3 \
                 --mmseqs2_module search \
                 --mmseqs2_ref "${longFnaX8}" \
                 --mmseqs2_qry "${longFnaX1}"
        ;;
    "prinseq")
        nextflow run ${nfDir}/lv1/prinseq.nf ${args} \
                 --prinseq_reads "${shortFnqGzPair}"
        ;;
    "spades")
        nextflow run ${nfDir}/lv1/spades.nf ${args} \
                 --spades_reads "${longFnqGzPair}"
        nextflow run ${nfDir}/lv1/spades.nf ${args} \
                 --spades_reads "${shortFnqGzPair}"
        ;;
    "virsorter")
        nextflow run ${nfDir}/lv1/virsorter.nf ${args} \
                 --virsorter_db "${virsorterDb}" \
                 --virsorter_mga "${virsorterMga}" \
                 --virsorter_db_type refseq \
                 --virsorter_aligner blast \
                 --virsorter_read "${longFnaX1}" \
        nextflow run ${nfDir}/lv1/virsorter.nf ${args} \
                 --virsorter_db "${virsorterDb}" \
                 --virsorter_mga "${virsorterMga}" \
                 --virsorter_db_type virome \
                 --virsorter_aligner diamond \
                 --virsorter_read "${longFnaX1}" \
        ;;
    "virsorter2")
        nextflow run ${nfDir}/lv1/virsorter2.nf ${args} \
                 --virsorter2_db "${virsorter2Db}" \
                 --virsorter2_read "${longFnaX1}"
        ;;
    "helloruby")
        nextflow run ${nfDir}/toys/helloruby.nf ${args} \
                 --helloruby_reads "${shortFnqGzPair}"
        ;;
    "*")
esac

if [ -f trace.${date}.txt ] ; then
    nfWorkDir="nfwork"
    workDir="${nfWorkDir}/work"

    while read line
    do
        hash=$(echo $line | awk '{print($2)}')
        if [ "${hash}" = "hash" ] ; then
            echo "hostname" > _
            continue
        fi
        hostname=$(head -n 1 ${workDir}/${hash}*/.command.log | awk '{print($5)}')
        echo ${hostname} >> _
    done<trace.${date}.txt

    paste _ trace.${date}.txt > __
    mv -f __ trace.${date}.txt
    rm -f _

    mv -f trace.${date}.txt trace_${test}.${date}.txt
fi

if [ -f report.${date}.html ] ; then
    mv -f report.${date}.html report_${test}.${date}.html
fi
    
if [ -f timeline.${date}.html ] ; then
    mv -f timeline.${date}.html timeline_${test}.${date}.html
fi

if [ -f dag.${date}.png ] ; then
    mv -f dag.${date}.png dag_${test}.${date}.png
fi

#rm -rf /dev/shm/${USER}
