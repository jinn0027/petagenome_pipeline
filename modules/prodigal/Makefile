MOD=prodigal

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prodigal -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prodigal -v



