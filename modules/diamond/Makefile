MOD=diamond

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} diamond help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} diamond version

