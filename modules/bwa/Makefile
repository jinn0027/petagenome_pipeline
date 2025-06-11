MOD=bwa

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bwa mem 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bwa 2>&1 | head -n 8


