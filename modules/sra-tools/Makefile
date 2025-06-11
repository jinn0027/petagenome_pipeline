MOD=sra-tools

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} fastq-dump -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} fastq-dump -h | tail -n 2

