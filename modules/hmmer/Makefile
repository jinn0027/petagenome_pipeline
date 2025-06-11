MOD=hmmer

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} jackhmmer -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} jackhmmer -h | head -n 2

