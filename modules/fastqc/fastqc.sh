#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/singularity exec ${MY_DIR}/fastqc.sif /opt/FastQC/fastqc --version


