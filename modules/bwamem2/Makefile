MOD=bwamem2

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bwa-mem2 mem 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bwa-mem2 2>&1 | head -n 8


