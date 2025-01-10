#!/bin/bash

date=$(date +%Y%m%d)

pushd ../external

git clone https://github.com/s-andrews/FastQC -b v0.12.1 --recursive
git clone https://github.com/marcelm/cutadapt -b v5.0 --recursive
git clone https://github.com/uwb-linux/prinseq --recursive
git clone https://github.com/ablab/spades -b v4.0.0 --recursive
git clone https://github.com/weizhongli/cdhit -b v4.8.1 --recursive
git clone https://github.com/biobakery/MetaPhlAn -b 4.1.1 --recursive
git clone https://github.com/simroux/VirSorter -b v1.0.6 --recursive
git clone https://github.com/EddyRivasLab/hmmer -b hmmer-3.4 --recursive
git clone https://github.com/arq5x/bedtools2 -b v2.31.1 --recursive
git clone https://github.com/shenwei356/seqkit -b v2.9.0 --recursive
git clone https://github.com/samtools/samtools -b 1.21 --recursive
git clone https://github.com/samtools/htslib -b 1.21 --recursive
git clone https://github.com/samtools/bcftools -b 1.21 --recursive
git clone https://github.com/BioInfoTools/BBMap -b v36.20 --recursive
git clone https://github.com/hyattpd/Prodigal -b v2.6.3 --recursive
git clone https://github.com/algbioi/ppsplus --recursive
git clone https://github.com/falcosecurity/falco -b v0.2.0 --recursive

for i in $(ls)
do
    if [ -d $i ] ; then
	echo $i
	tar cvfz ${i}.${date}.tar.gz ${i}
	rm -rf $i
    fi
done

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
mv ncbi-blast-2.16.0+-x64-linux.tar.gz ncbi-blast-2.16.0+-x64-linux.${date}.tar.gz

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
mv fastqc_v0.12.1.zip fastqc_v0.12.1.${date}.zip

popd
