#!/bin/bash

mkdir -p ~/bin

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
