#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

/usr/local/bin/apptainer exec --fakeroot virsorter.sif wrapper_phage_contigs_sorter_iPlant.pl -h







