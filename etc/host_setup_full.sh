#!/bin/bash

sudo dnf -y groupinstall 'Development Tools'
sudo dnf -y install \
     openssl-devel \
     libuuid-devel \
     libseccomp-devel \
     wget \
     git \
     squashfs-tools \
     glib2-devel \
     fuse3-devel \
     fakeroot \
     cryptsetup \
     squashfuse \
     libsubid-dev \
     pigz \
     bzip2 \
     pbzip2

#sudo dnf --enablerepo=devel install -y \
#     shadow-utils-subid-devel

sudo dnf install -y \
     shadow-utils-subid

# apptainer/singularityにgo必要なのでインストール
if ! (type "go" > /dev/null 2>&1); then
    if [ ! -d /usr/local/go ] ; then
        if [ ! -f go1.23.3.linux-amd64.tar.gz ] ; then
            wget https://go.dev/dl/go1.23.3.linux-amd64.tar.gz
        fi
        sudo tar -I pigz -C /usr/local -xf go1.23.3.linux-amd64.tar.gz
    fi
    export PATH=/usr/local/go/bin:$PATH
fi

# apptainerインストール
if ! (type "apptainer" > /dev/null 2>&1); then
    if [ ! -d apptainer ] ; then
        git clone https://github.com/apptainer/apptainer.git
    fi
    cd apptainer
    ./mconfig
    make -C ./builddir
    sudo make -C ./builddir install
    cd ..
fi

# singularityインストール
#if ! (type "singularity" > /dev/null 2>&1); then
#    if [ ! -d singularity ] ; then
#        git clone https://github.com/sylabs/singularity.git -b v4.2.2 --recursive
#    fi
#    cd singularity
#    ./mconfig
#    make -C ./builddir
#    sudo make -C ./builddir install
#    cd ..
#fi

# nextflowにjava必要なのでインストール
if ! (type "java" > /dev/null 2>&1); then
    curl -s https://get.sdkman.io | bash
    source ~/.sdkman/bin/sdkman-init.sh
    sdk install java 17.0.10-tem
    java --version
    #openjdk 17.0.10 2024-01-16
    #OpenJDK Runtime Environment Temurin-17.0.10+7 (build 17.0.10+7)
    #OpenJDK 64-Bit Server VM Temurin-17.0.10+7 (build 17.0.10+7, mixed mode, sharing)
fi

# nextflowインストール
if ! (type "nextflow" > /dev/null 2>&1); then
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin
    nextflow -v
    #nextflow version 24.10.2.5932
fi

###############################################

# SElinuxをpermissiveに設定。
# これをやらないとapptainerでfakerootでdnfをやろうとするとエラーが発生する。
sudo setenforce 0

# apptainerのキャッシュ先を/dev/shmに
mkdir -p /dev/shm/.apptainer
export APPTAINER_CACHEDIR=/dev/shm/.apptainer

