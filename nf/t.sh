#!/bin/bash

nextflow run main.nf --reads "../modules/test/s_6_{1,2}.fastq.gz" --threads 4 --petagenomeDir=$(pwd)/..

#nextflow run fastp.nf --reads "../modules/test/s_6_{1,2}.fastq.gz" --threads $(nproc) --petagenomeDir=$(pwd)/..

