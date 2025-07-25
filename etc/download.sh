#!/bin/bash

#date=$(date +%Y%m%d)

pushd ../external

#if [ ! -f FastQC.tar.gz ] ; then
#    git clone https://github.com/s-andrews/FastQC -b v0.12.1 --recursive
#fi

if [ ! -f cutadapt.tar.gz ] ; then
    git clone https://github.com/marcelm/cutadapt -b v5.0 --recursive
fi

if [ ! -f prinseq.tar.gz ] ; then
    git clone https://github.com/uwb-linux/prinseq --recursive
fi

if [ ! -f spades.tar.gz ] ; then
    git clone https://github.com/ablab/spades -b v4.0.0 --recursive
fi

if [ ! -f cdhit.tar.gz ] ; then
    git clone https://github.com/weizhongli/cdhit -b V4.8.1 --recursive
fi

if [ ! -f MetaPhlAn.tar.gz ] ; then
    git clone https://github.com/biobakery/MetaPhlAn -b 4.1.1 --recursive
fi

if [ ! -f VirSorter.tar.gz ] ; then
    git clone https://github.com/simroux/VirSorter -b v1.0.6 --recursive
fi

if [ ! -f VirSorter2.tar.gz ] ; then
    git clone https://github.com/jiarong/VirSorter2 -b v2.2.4 --recursive
fi

if [ ! -f hmmer.tar.gz ] ; then
    git clone https://github.com/EddyRivasLab/hmmer -b infernal-1.1.5 --recursive
    pushd hmmer
    git clone https://github.com/EddyRivasLab/easel -b infernal-1.1.5 --recursive
    popd
fi

if [ ! -f bedtools2.tar.gz ] ; then
    git clone https://github.com/arq5x/bedtools2 -b v2.31.1 --recursive
fi

if [ ! -f seqkit.tar.gz ] ; then
    git clone https://github.com/shenwei356/seqkit -b v2.9.0 --recursive
fi

if [ ! -f samtools.tar.gz ] ; then
    git clone https://github.com/samtools/samtools -b 1.21 --recursive
fi

if [ ! -f htslib.tar.gz ] ; then
    git clone https://github.com/samtools/htslib -b 1.21 --recursive
fi

if [ ! -f bcftools.tar.gz ] ; then
    git clone https://github.com/samtools/bcftools -b 1.21 --recursive
fi

if [ ! -f BBMap.tar.gz ] ; then
    git clone https://github.com/BioInfoTools/BBMap -b v36.20 --recursive
fi

if [ ! -f Prodigal.tar.gz ] ; then
    git clone https://github.com/hyattpd/Prodigal -b v2.6.3 --recursive
fi

if [ ! -f ppsplus.tar.gz ] ; then
    git clone https://github.com/algbioi/ppsplus --recursive
fi

if [ ! -f falco.tar.gz ] ; then
    git clone https://github.com/smithlabcode/falco -b v1.2.5 --recursive
fi

if [ ! -f bwa.tar.gz ] ; then
    git clone https://github.com/lh3/bwa -b v0.7.18 --recursive
fi

if [ ! -f bwa-mem2.tar.gz ] ; then
    git clone https://github.com/bwa-mem2/bwa-mem2 -b v2.2.1 --recursive
fi

if [ ! -f diamond.0.9.36.tar.gz ] ; then
    git clone https://github.com/bbuchfink/diamond -b v0.9.36 --recursive
    mv diamond diamond.0.9.36
fi

if [ ! -f diamond.tar.gz ] ; then
    git clone https://github.com/bbuchfink/diamond -b v2.1.10 --recursive
fi

if [ ! -f bowtie.tar.gz ] ; then
    git clone https://github.com/BenLangmead/bowtie -b v1.3.1 --recursive
fi

if [ ! -f bowtie2.tar.gz ] ; then
    git clone https://github.com/BenLangmead/bowtie2 -b v2.5.4 --recursive
fi

if [ ! -f sra-tools.tar.gz ] ; then
    git clone https://github.com/ncbi/sra-tools -b 3.2.0 --recursive
fi

if [ ! -f ncbi-vdb.tar.gz ] ; then
    git clone https://github.com/ncbi/ncbi-vdb -b 3.2.0 --recursive
fi

if [ ! -f megahit.tar.gz ] ; then
    git clone https://github.com/voutcn/megahit -b v1.2.9 --recursive
fi

if [ ! -f minimap2.tar.gz ] ; then
    git clone https://github.com/lh3/minimap2 -b v2.28 --recursive
fi

if [ ! -f MMseqs2.tar.gz ] ; then
    git clone https://github.com/soedinglab/MMseqs2.git -b 17-b804f --recursive
fi

#if [ ! -f yasm.tar.gz ] ; then
#    git clone https://github.com/yasm/yasm -b v1.3.0 --recursive
#fi

#if [ ! -f nasm.tar.gz ] ; then
#    git clone https://github.com/netwide-assembler/nasm -b nasm-2.16.03 --recursive
#fi

#if [ ! -f isa-l.tar.gz ] ; then
#    git clone https://github.com/intel/isa-l -b v2.31.1 --recursive
#fi

#if [ ! -f libdeflate.tar.gz ] ; then
#    git clone https://github.com/ebiggers/libdeflate -b v1.23 --recursive
#fi

if [ ! -f fastp.tar.gz ] ; then
    git clone https://github.com/OpenGene/fastp -b v0.24.0 --recursive
fi

for i in $(ls)
do
    if [ -d $i ] ; then
	if [ "$i" == "metaphlan_db" ] ; then
	    continue
	fi
	if [ "$i" == "virsorter-data" ] ; then
	    continue
	fi
	if [ "$i" == "virsorter2-data" ] ; then
	    continue
	fi
        echo $i
        tar -I pigz -cvf ${i}.tar.gz ${i}
        rm -rf $i
    fi
done

if [ ! -f seqkit_linux_amd64.2.9.0.tar.gz ] ; then
    wget -O seqkit_linux_amd64.2.9.0.tar.gz  https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz
fi

# MiniConda setup script
if [ ! -f Miniconda3-py39_24.9.2-0-Linux-x86_64.sh ] ; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.9.2-0-Linux-x86_64.sh
fi

# Blast+ binary
if [ ! -f ncbi-blast-2.16.0+-x64-linux.tar.gz ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
fi

# Blast+ source
if [ ! -f ncbi-blast-2.16.0+-src.tar.gz ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-src.tar.gz
fi

# FastQC binary
if [ ! -f fastqc_v0.12.1.zip ] ; then
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
fi

# MetaPhlAn database @ 2025/1/29

if [ ! -d metaphlan_db ] ; then
    metaphlan_index=mpa_vJun23_CHOCOPhlAnSGB_202403
    echo ${metaphlan_index} > mpa_latest

    if [ ! -f ${metaphlan_index}.md5 ] ; then
        wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${metaphlan_index}.md5
    fi

    if [ ! -f ${metaphlan_index}.tar.gz ] ; then
        wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${metaphlan_index}.tar
        gzip ${metaphlan_index}.tar
    fi

    if [ ! -f ${metaphlan_index}_bt2.md5 ] ; then
        wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/${metaphlan_index}_bt2.md5
    fi

    if [ ! -f ${metaphlan_index}_bt2.tar.gz ] ; then
        wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/${metaphlan_index}_bt2.tar
        gzip ${metaphlan_index}_bt2.tar
    fi

    mkdir -p metaphlan_db
    cp mpa_latest metaphlan_db
    tar -I pigz -xvf ${metaphlan_index}.tar.gz -C metaphlan_db --no-same-owner
    tar -I pigz -xvf ${metaphlan_index}_bt2.tar.gz -C metaphlan_db --no-same-owner
    pushd metaphlan_db
        pbunzip2 *.bz2
    popd
    chmod -R 777 metaphlan_db
fi

# VirSorter Metagenome Annotator @ 2012/10/7 (retrieved at 2025/1/22)
if [ ! -f mga_x86_64.tar.gz ] ; then
    wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
fi

# VirSorter database @ 2025/1/16
if [ ! -f virsorter-data-v2.tar.gz ] ; then
    wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
fi

# VirSorter2 database @ 2025/1/22
if [ ! -d virsorter2-data ] ; then
    if [ ! -f db.tar.gz ] ; then
        wget https://osf.io/v46sc/download -O db.tar.gz
    fi

    if [ ! -d virsorter2-data ] ; then
        tar -I pigz -xvf db.tar.gz
        mv db virsorter2-data
    fi
fi

# VirSorter2 test @ 2025/1/22
if [ ! -f 8seq.fa ] ; then
    wget https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
fi

# BBMap v39.15
if [ ! -f BBMap_39.15.tar.gz ] ; then
    wget https://sourceforge.net/projects/bbmap/files/BBMap_39.15.tar.gz/download -O BBMap_39.15.tar.gz
fi

# Eclipse jar for compiling BBMap
#if [ ! -f ecj-4.5.2.jar ] ; then
#    wget http://www.eclipse.org/downloads/download.php?file=/eclipse/downloads/drops4/R-4.5.2-201602121500/ecj-4.5.2.jar -O ecj-4.5.2.jar
#fi

# hmmer Pfam database vir interpro @ 2025/1/31 via interpro( https://www.ebi.ac.uk/interpro/download/pfam/ )
# ref) https://qiita.com/116ryusei/items/5fa0f1d8291c046cffe7
if [ ! -f Pfam-A.hmm.gz ] ; then
    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
fi

# modify BBMap
if [ ! -f BBMap.modified.tar.gz ] ; then
    tar -I pigz -xvf BBMap.tar.gz
    sed -i -e 's#/usr/common/usg/hpc/openmpi/gnu4.6/sge/1.8.1/ib_2.1-1.0.0/lib/mpi.jar#/usr/lib64/openmpi/lib/mpi.jar#g' \
           -e 's#compiler="${jcompiler}"##g' BBMap/build.xml
    tar -I pigz -cvf BBMap.modified.tar.gz BBMap
    rm -rf BBMap
fi

# modify BBMap_39.15
if [ ! -f BBMap_39.15.modified.tar.gz ] ; then
    tar -I pigz -xvf BBMap_39.15.tar.gz
    sed -i -e 's#/usr/common/usg/hpc/openmpi/gnu4.6/sge/1.8.1/ib_2.1-1.0.0/lib/mpi.jar#/usr/lib64/openmpi/lib/mpi.jar#g' \
           -e 's#compiler="${jcompiler}"##g' bbmap/build.xml
    tar -I pigz -cvf BBMap_39.15.modified.tar.gz bbmap
    rm -rf bbmap
fi

# update virsorter database
# see modules/virsorter/memo
if [ ! -f virsorter-data-v2.updated.tar.gz ] ; then
    pushd ../modules/common
    make gcc.sif
    popd
    pushd ../modules/virsorter
    make virsorter_update_db.sbx
    make db_updated
    popd
fi

if [ ! -d virsorter-data ] ; then
    tar -I pigz -xvf virsorter-data-v2.updated.tar.gz
fi

if [ ! -f mga_linux_ia64 ] ; then
    tar xvfz mga_x86_64.tar.gz mga_linux_ia64
fi

popd # external
