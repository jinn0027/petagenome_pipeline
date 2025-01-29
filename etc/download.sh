#!/bin/bash

#date=$(date +%Y%m%d)

pushd ../external

#git clone https://github.com/s-andrews/FastQC -b v0.12.1 --recursive
git clone https://github.com/marcelm/cutadapt -b v5.0 --recursive
git clone https://github.com/uwb-linux/prinseq --recursive
git clone https://github.com/ablab/spades -b v4.0.0 --recursive
git clone https://github.com/weizhongli/cdhit -b V4.8.1 --recursive
git clone https://github.com/biobakery/MetaPhlAn -b 4.1.1 --recursive
git clone https://github.com/biobakery/MetaPhlAn2 --recursive
git clone https://github.com/simroux/VirSorter -b v1.0.6 --recursive
git clone https://github.com/jiarong/VirSorter2 -b v2.2.4 --recursive
#git clone https://github.com/EddyRivasLab/hmmer -b hmmer-3.4 --recursive
git clone https://github.com/EddyRivasLab/hmmer -b infernal-1.1.5 --recursive
pushd hmmer
git clone https://github.com/EddyRivasLab/easel -b infernal-1.1.5 --recursive
popd
git clone https://github.com/arq5x/bedtools2 -b v2.31.1 --recursive
git clone https://github.com/shenwei356/seqkit -b v2.9.0 --recursive
wget -O seqkit_linux_amd64.2.9.0.tar.gz  https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz

git clone https://github.com/samtools/samtools -b 1.21 --recursive
git clone https://github.com/samtools/htslib -b 1.21 --recursive
git clone https://github.com/samtools/bcftools -b 1.21 --recursive

git clone https://github.com/BioInfoTools/BBMap -b v36.20 --recursive
#git clone https://github.com/BioInfoTools/BBMap.git --recursive
#pushd BBMap
#  git checkout 70b24bd416b158e1738f8b70682f3bd2e656ed60
#popd

git clone https://github.com/hyattpd/Prodigal -b v2.6.3 --recursive
git clone https://github.com/algbioi/ppsplus --recursive
git clone https://github.com/smithlabcode/falco -b v1.2.5 --recursive
git clone https://github.com/lh3/bwa -b v0.7.18 --recursive

# diamond with blast seems not possible on v2.1.10
git clone https://github.com/bbuchfink/diamond -b v2.1.10 --recursive
# 0.9.14 has an installation bug : http://diamondsearch.org/forums/index.php?threads/diamond-installation-error.4/
# 0.9.36 can be installed successfully but yet has a version mismatch error of database for virsorter.
#  as : Error: Database was built with an older version of Diamond and is incompatible.
#git clone https://github.com/bbuchfink/diamond -b v0.9.14 --recursive diamond.0.9.14
#git clone https://github.com/bbuchfink/diamond -b v0.9.36 --recursive diamond.0.9.36
git clone https://github.com/bbuchfink/diamond -b v0.9.30 --recursive diamond.0.9.25

git clone https://github.com/BenLangmead/bowtie -b v1.3.1 --recursive
git clone https://github.com/BenLangmead/bowtie2 -b v2.5.4 --recursive
git clone https://github.com/ncbi/sra-tools -b 3.2.0 --recursive
git clone https://github.com/ncbi/ncbi-vdb -b 3.2.0 --recursive
git clone https://github.com/voutcn/megahit -b v1.2.9 --recursive

#git clone https://github.com/ncbi/ncbi-cxx-toolkit-public.git -b release-28.0.12 --recursive

for i in $(ls)
do
    if [ -d $i ] ; then
	echo $i
	tar cvfz ${i}.tar.gz ${i}
	rm -rf $i
    fi
done

# MiniConda setup script
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.9.2-0-Linux-x86_64.sh

# Blast+ binary
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
# Blast+ source
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.16.0+-src.tar.gz

# FastQC binary
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip

# MetaPhlAn2 database @ 2025/1/22
wget https://figshare.com/ndownloader/articles/6200807/versions/1 -O mpa_v20_m200.zip

# VirSorter Metagenome Annotator @ 2012/10/7 (retrieved at 2025/1/22)
wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
# VirSorter database @ 2025/1/16
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
# virsorter-data-v2.updated.tar.gz : VirSorterr database updated
# see modules/virsorter/memo

# VirSorter2 database @ 2025/1/22
wget https://osf.io/v46sc/download -O db.tar.gz
# VirSorter2 test @ 2025/1/22
wget https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa


# Eclipse jar for compiling BBMap
wget http://www.eclipse.org/downloads/download.php?file=/eclipse/downloads/drops4/R-4.5.2-201602121500/ecj-4.5.2.jar -O ecj-4.5.2.jar

popd
