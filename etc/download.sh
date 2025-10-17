#!/bin/bash

#date=$(date +%Y%m%d)

NEED_MODULES=$(pwd)/need_modules.txt
DIR_MODULES=$(pwd)/../modules
DIR_EXTERNAL=$(pwd)/../external

function check_module () {
    module=$1
    if [ ! -f ${NEED_MODULES} ] || grep -w -q ${module} ${NEED_MODULES} ; then
        return 0
    else
        return 1
    fi
}

function need_make () {
    file=$1
    if [ -f ${DIR_EXTERNAL}/${file} ] ; then
	return 1
    fi
    if [ ! -f ${NEED_MODULES} ] ; then
        return 0
    fi
    while read line ; do
        if [ "${line}" = "" ] ; then
            continue
        fi
        dir=${DIR_MODULES}/${line}
        if [ -d ${dir} ] ; then
            for def in $(ls ${dir}/*.def) ; do
		if grep -w -q ${file} ${def} ; then
                    return 0
		fi
            done
        fi
    done<${NEED_MODULES}
    return 1
}

pushd ../external

############ MiniConda setup script

if [ ! -f Miniconda3-py39_24.9.2-0-Linux-x86_64.sh ] ; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.9.2-0-Linux-x86_64.sh
fi

############ source code or binary code

if need_make BBMap.tar.gz || need_make BBMap.modified.tar.gz ; then
    git clone https://github.com/BioInfoTools/BBMap -b v36.20 --recursive
    tar -I pigz -cvf BBMap.tar.gz BBMap
    rm -rf BBMap
fi

if need_make BBMap_39.15.tar.gz || need_make BBMap_39.15.modified.tar.gz ; then
    wget https://sourceforge.net/projects/bbmap/files/BBMap_39.15.tar.gz/download -O BBMap_39.15.tar.gz
fi

# modify BBMap
if need_make BBMap.modified.tar.gz ; then
    tar -I pigz -xvf BBMap.tar.gz
    sed -i -e 's#/usr/common/usg/hpc/openmpi/gnu4.6/sge/1.8.1/ib_2.1-1.0.0/lib/mpi.jar#/usr/lib64/openmpi/lib/mpi.jar#g' \
           -e 's#compiler="${jcompiler}"##g' BBMap/build.xml
    tar -I pigz -cvf BBMap.modified.tar.gz BBMap
    rm -rf BBMap
fi

# modify BBMap_39.15
if need_make BBMap_39.15.modified.tar.gz ; then
    tar -I pigz -xvf BBMap_39.15.tar.gz
    sed -i -e 's#/usr/common/usg/hpc/openmpi/gnu4.6/sge/1.8.1/ib_2.1-1.0.0/lib/mpi.jar#/usr/lib64/openmpi/lib/mpi.jar#g' \
           -e 's#compiler="${jcompiler}"##g' bbmap/build.xml
    tar -I pigz -cvf BBMap_39.15.modified.tar.gz bbmap
    rm -rf bbmap
fi

if need_make bcftools.tar.gz ; then
    git clone https://github.com/samtools/bcftools -b 1.21 --recursive
    tar -I pigz -cvf bcftools.tar.gz bcftools
    rm -rf bcftools
fi

if need_make bedtools2.tar.gz ; then
    git clone https://github.com/arq5x/bedtools2 -b v2.31.1 --recursive
    tar -I pigz -cvf bedtools2.tar.gz bedtools2
    rm -rf bedtools2
fi

if need_make bwa.tar.gz ; then
    git clone https://github.com/lh3/bwa -b v0.7.18 --recursive
    tar -I pigz -cvf bwa.tar.gz bwa
    rm -rf bwa
fi

if need_make bwa-mem2.tar.gz ; then
    git clone https://github.com/bwa-mem2/bwa-mem2 -b v2.2.1 --recursive
    tar -I pigz -cvf bwa-mem2.tar.gz bwa-mem2
    rm -rf bwa-mem2
fi

if need_make bowtie.tar.gz ; then
    git clone https://github.com/BenLangmead/bowtie -b v1.3.1 --recursive
    tar -I pigz -cvf bowtie.tar.gz bowtie
    rm -rf biwtie
fi

if need_make bowtie2.tar.gz ; then
    git clone https://github.com/BenLangmead/bowtie2 -b v2.5.4 --recursive
    tar -I pigz -cvf bowtie2.tar.gz bowtie2
    rm -rf biwtie2
fi

if need_make cdhit.tar.gz ; then
    git clone https://github.com/weizhongli/cdhit -b V4.8.1 --recursive
    tar -I pigz -cvf cdhit.tar.gz cdhit
    rm -rf cdhit
fi

if need_make cutadapt.tar.gz ; then
    git clone https://github.com/marcelm/cutadapt -b v5.0 --recursive
    tar -I pigz -cvf cutadapt.tar.gz cutadapt
    rm -rf cutadapt
fi

if need_make diamond.tar.gz ; then
    git clone https://github.com/bbuchfink/diamond -b v2.1.10 --recursive
    tar -I pigz -cvf diamond.tar.gz diamond
    rm -rf diamond
fi

if need_make diamond.0.9.36.tar.gz ; then
    git clone https://github.com/bbuchfink/diamond -b v0.9.36 --recursive
    mv diamond diamond.0.9.36
    tar -I pigz -cvf diamond.0.9.36.tar.gz diamond.0.9.36
    rm -rf diamond.0.9.36
fi

if need_make falco.tar.gz ; then
    git clone https://github.com/smithlabcode/falco -b v1.2.5 --recursive
    tar -I pigz -cvf falco.tar.gz falco
    rm -rf falco
fi

if need_make fastp.tar.gz ; then
    git clone https://github.com/OpenGene/fastp -b v0.24.0 --recursive
    tar -I pigz -cvf fastp.tar.gz fastp
    rm -rf fastp
fi

if need_make fastqc_v0.12.1.zip ; then
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
fi

if need_make hmmer.tar.gz ; then
    git clone https://github.com/EddyRivasLab/hmmer -b infernal-1.1.5 --recursive
    pushd hmmer
        git clone https://github.com/EddyRivasLab/easel -b infernal-1.1.5 --recursive
    popd
    tar -I pigz -cvf hmmer.tar.gz hmmer
    rm -rf hmmer
fi

if need_make htslib.tar.gz ; then
    git clone https://github.com/samtools/htslib -b 1.21 --recursive
    tar -I pigz -cvf htslib.tar.gz htslib
    rm -rf htslib
fi

if need_make jloda.tar.gz ; then
    git clone https://github.com/husonlab/jloda --recursive
    pushd jloda
        git checkout 27ec9d8
    popd
    tar -I pigz -cvf jloda.tar.gz jloda
    rm -rf jloda
fi

if need_make MetaCarvel.tar.gz ; then
    git clone https://github.com/marbl/MetaCarvel --recursive
    pushd MetaCarvel
        git checkout d4b5408
    popd
    tar -I pigz -cvf MetaCarvel.tar.gz MetaCarvel
    rm -rf MetaCarvel
fi

if need_make megahit.tar.gz ; then
    git clone https://github.com/voutcn/megahit -b v1.2.9 --recursive
    tar -I pigz -cvf megahit.tar.gz megahit
    rm -rf megahit
fi

if need_make megan-ce.tar.gz ; then
    git clone https://github.com/husonlab/megan-ce --recursive
    pushd megan-ce
        git checkout 9b36523
    popd
    tar -I pigz -cvf megan-ce.tar.gz megan-ce
    rm -rf megan-ce
fi

if need_make MEGAN_Community_unix_6_25_10.sh ; then
    wget https://software-ab.cs.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_25_10.sh
fi

if need_make MetaPhlAn.tar.gz ; then
    git clone https://github.com/biobakery/MetaPhlAn -b 4.1.1 --recursive
    tar -I pigz -cvf MetaPhlAn.tar.gz MetaPhlAn
    rm -rf MetaPhlAn
fi

if need_make minimap2.tar.gz ; then
    git clone https://github.com/lh3/minimap2 -b v2.28 --recursive
    tar -I pigz -cvf minimap2.tar.gz minimap2
    rm -rf minimap2
fi

if need_make MMseqs2.tar.gz ; then
    git clone https://github.com/soedinglab/MMseqs2.git -b 17-b804f --recursive
    tar -I pigz -cvf MMseqs2.tar.gz MMseqs2
    rm -rf MMseqs2
fi

if need_make ncbi-blast-2.16.0+-src.tar.gz ; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-src.tar.gz
fi

if need_make ncbi-blast-2.16.0+-x64-linux.tar.gz ; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
fi

if need_make openjfx-24.0.2_linux-aarch64_bin-sdk.zip ; then
    wget https://download2.gluonhq.com/openjfx/24.0.2/openjfx-24.0.2_linux-aarch64_bin-sdk.zip
fi

if need_make ncbi-vdb.tar.gz ; then
    git clone https://github.com/ncbi/ncbi-vdb -b 3.2.0 --recursive
    tar -I pigz -cvf ncbi-vdb.tar.gz ncbi-vdb
    rm -rf ncbi-vdb
fi

if need_make prinseq.tar.gz ; then
    git clone https://github.com/uwb-linux/prinseq --recursive
    tar -I pigz -cvf prinseq.tar.gz prinseq
    rm -rf prinseq
fi

if need_make Prodigal.tar.gz ; then
    git clone https://github.com/hyattpd/Prodigal -b v2.6.3 --recursive
    tar -I pigz -cvf Prodigal.tar.gz Prodigal
    rm -rf Prodigal
fi

if need_make samtools.tar.gz ; then
    git clone https://github.com/samtools/samtools -b 1.21 --recursive
    tar -I pigz -cvf samtools.tar.gz samtools
    rm -rf samtools
fi

if need_make seqkit.tar.gz ; then
    git clone https://github.com/shenwei356/seqkit -b v2.9.0 --recursive
    tar -I pigz -cvf seqkit.tar.gz seqkit
    rm -rf seqkit
fi

if need_make seqkit_linux_amd64.2.9.0.tar.gz ; then
    wget -O seqkit_linux_amd64.2.9.0.tar.gz  https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz
fi

if need_make SOAPdenovo2.tar.gz ; then
    git clone  https://github.com/aquaskyline/SOAPdenovo2 --recursive
    pushd SOAPdenovo2
        git checkout 7560c0b
    popd
    # 以下はコンパイル通らない
    #git clone  https://github.com/aquaskyline/SOAPdenovo2 -b r242 --recursive
    tar -I pigz -cvf SOAPdenovo2.tar.gz SOAPdenovo2
    rm -rf SOAPdenovo2
fi

if need_make spades.tar.gz ; then
    git clone https://github.com/ablab/spades -b v4.0.0 --recursive
    tar -I pigz -cvf spades.tar.gz spades
    rm -rf spades
fi

if need_make sra-tools.tar.gz ; then
    git clone https://github.com/ncbi/sra-tools -b 3.2.0 --recursive
    tar -I pigz -cvf sra-tools.tar.gz sra-tools
    rm -rf sra-tools
fi

if need_make sspace_basic.tar.gz ; then
    git clone https://github.com/nsoranzo/sspace_basic.git -b v2.1.1 --recursive
    tar -I pigz -cvf sspace_basic.tar.gz sspace_basic
    rm -rf sspace_basic
fi

if need_make VirSorter.tar.gz ; then
    git clone https://github.com/simroux/VirSorter -b v1.0.6 --recursive
    tar -I pigz -cvf VirSorter.tar.gz VirSorter
    rm -rf VirSorter
fi

if need_make VirSorter2.tar.gz ; then
    git clone https://github.com/jiarong/VirSorter2 -b v2.2.4 --recursive
    tar -I pigz -cvf VirSorter2.tar.gz VirSorter2
    rm -rf VirSorter2
fi

############ database

# hmmer Pfam database vir interpro @ 2025/1/31 via interpro( https://www.ebi.ac.uk/interpro/download/pfam/ )
# ref) https://qiita.com/116ryusei/items/5fa0f1d8291c046cffe7
if check_module hmmer && [ ! -f Pfam-A.hmm.gz ] ; then
    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
fi

# MetaPhlAn database @ 2025/1/29
if check_module metaphlan && [ ! -d metaphlan_db ] ; then
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
if check_module virsorter && [ ! -f mga_x86_64.tar.gz ] ; then
    wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
fi

# VirSorter database @ 2025/1/16
if check_module virsorter && [ ! -f virsorter-data-v2.tar.gz ] ; then
    wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
fi

# update VirSorter database
# see modules/virsorter/memo
if check_module virsorter && [ ! -f virsorter-data-v2.updated.tar.gz ] ; then
    pushd ../modules/common
        make gcc.sif
    popd
    pushd ../modules/virsorter
        make virsorter_update_db.sbx
        make db_updated
    popd
fi

# extract updated VirSorter database
if check_module virsorter && [ ! -d virsorter-data ] ; then
    tar -I pigz -xvf virsorter-data-v2.updated.tar.gz
fi

# VirSorter2 database @ 2025/1/22
if check_module virsorter2 && [ ! -d virsorter2-data ] ; then
    if [ ! -f db.tar.gz ] ; then
        wget https://osf.io/v46sc/download -O db.tar.gz
    fi

    if [ ! -d virsorter2-data ] ; then
        tar -I pigz -xvf db.tar.gz
        mv db virsorter2-data
    fi
fi

popd # external
