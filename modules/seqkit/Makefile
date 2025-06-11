MOD=seqkit

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} seqkit --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} seqkit version

