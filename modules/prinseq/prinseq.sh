#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/apptainer exec ${MY_DIR}/prinseq.sif perl /opt/prinseq/prinseq-lite.pl -h



