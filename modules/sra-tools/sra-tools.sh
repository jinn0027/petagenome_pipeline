#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/apptainer exec sra-tools.sif fastq-dump -h





