#!/bin/bash

# singularityにgo必要なのでインストール
wget https://go.dev/dl/go1.23.3.linux-amd64.tar.gz
sudo tar -C /usr/local -xzf go1.23.3.linux-amd64.tar.gz
export PATH=/usr/local/go/bin:$PATH

sudo yum -y groupinstall 'Development Tools'
sudo yum -y install openssl-devel libuuid-devel libseccomp-devel \
                    wget squashfs-tools glib2-devel-2.68.4-14.el9_4.1.x86_64 \
                    fuse3-devel-3.10.2-9.el9.x86_64
wget https://github.com/sylabs/singularity/releases/download/v4.2.2/singularity-ce-4.2.2.tar.gz

# singularityインストール

tar xvfz singularity-ce-4.2.2.tar.gz
cd singularity-ce-4.2.2
./mconfig
make -C ./builddir
sudo make -C ./builddir install
cd ..

# nextflowにjava必要なのでインストール

curl -s https://get.sdkman.io | bash
source ~/.sdkman/bin/sdkman-init.sh
sdk install java 17.0.10-tem
java --version
#openjdk 17.0.10 2024-01-16
#OpenJDK Runtime Environment Temurin-17.0.10+7 (build 17.0.10+7)
#OpenJDK 64-Bit Server VM Temurin-17.0.10+7 (build 17.0.10+7, mixed mode, sharing)

# nextflowインストール
curl -s https://get.nextflow.io | bash
chmod +x nextflow

sudo mv nextflow /usr/local/bin
nextflow -v
#nextflow version 24.10.2.5932
