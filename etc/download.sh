#!/bin/bash

#date=$(date +%Y%m%d)

pushd ../external

#git clone https://github.com/s-andrews/FastQC -b v0.12.1 --recursive
git clone https://github.com/marcelm/cutadapt -b v5.0 --recursive
git clone https://github.com/uwb-linux/prinseq --recursive
git clone https://github.com/ablab/spades -b v4.0.0 --recursive
git clone https://github.com/weizhongli/cdhit -b V4.8.1 --recursive
git clone https://github.com/biobakery/MetaPhlAn -b 4.1.1 --recursive
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

git clone https://github.com/hyattpd/Prodigal -b v2.6.3 --recursive
git clone https://github.com/algbioi/ppsplus --recursive
git clone https://github.com/smithlabcode/falco -b v1.2.5 --recursive
git clone https://github.com/lh3/bwa -b v0.7.18 --recursive
git clone https://github.com/bbuchfink/diamond -b v2.1.10 --recursive
git clone https://github.com/BenLangmead/bowtie -b v1.3.1 --recursive
git clone https://github.com/BenLangmead/bowtie2 -b v2.5.4 --recursive
git clone https://github.com/ncbi/sra-tools -b 3.2.0 --recursive
git clone https://github.com/ncbi/ncbi-vdb -b 3.2.0 --recursive

for i in $(ls)
do
    if [ -d $i ] ; then
	echo $i
	tar cvfz ${i}.tar.gz ${i}
	rm -rf $i
    fi
done

wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.9.2-0-Linux-x86_64.sh

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip

wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
#tar xvfz virsorter-data-v2.tar.gz

popd
