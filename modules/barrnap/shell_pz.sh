#!/bin/bash

MOD=barrnap

SIF=${MOD}.sif
OVL=${MOD}.ovl
SBX=${MOD}.sbx
WORK=$(pwd)

SINGULARITY=apptainer

${SINGULARITY} shell --no-home \
                     --bind ${WORK}:/opt/work \
                     --pwd /opt \
                     --shell /bin/bash \
                     ${SIF}

