#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/singularity exec ${MY_DIR}/prinseq.sif /opt/prinseq


