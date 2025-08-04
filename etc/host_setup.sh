#!/bin/bash

if groups | grep -q wheel ; then
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
        pigz \
        bzip2 \
        pbzip2

    #sudo dnf --enablerepo=devel install -y \
    #    shadow-utils-subid-devel

    sudo dnf install -y \
        shadow-utils-subid

    # SElinuxをpermissiveに設定。
    # これをやらないとapptainerでfakerootでdnfをやろうとするとエラーが発生する場合がある。
    sudo setenforce 0
fi

mkdir -p ~/bin

# nextflowにjava必要

need_java=1
if (type "java" > /dev/null 2>&1); then
    ver=$(java -version 2>&1 | grep "version" | awk -F '"' '{print $2}')
    major=$(echo ${ver} | awk -F '.' '{print($1)}')
    if [ ${major} -ge 17 ] ; then
	need_java=0
    fi
fi

if [ ${need_java} = "1" ] ; then
    curl -s https://get.sdkman.io | bash
    source ~/.sdkman/bin/sdkman-init.sh
    sdk install java 17.0.10-tem
    java --version
fi

if ! (type "apptainer" > /dev/null 2>&1); then
    curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | bash -s - ~
fi

if ! (type "nextflow" > /dev/null 2>&1); then
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow ~/bin
fi

if ! (type "apptainer" > /dev/null 2>&1) || ! (type "nextflow" > /dev/null 2>&1) ; then
    export PATH=~/bin:$PATH
fi
export PETAGENOME_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
