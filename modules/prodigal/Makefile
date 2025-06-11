MOD=prodigal

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prodigal -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prodigal -v



