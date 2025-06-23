#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/apptainer exec metaphlan2.sif metaphlan2.py --help




