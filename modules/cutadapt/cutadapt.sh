#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/singularity exec --fakeroot ${MY_DIR}/cutadapt.sif cutadapt -h



