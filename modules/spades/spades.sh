#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

apptainer exec ${MY_DIR}/spades.sif spades.py -h




