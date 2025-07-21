MOD=spades

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} spades.py -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} spades.py -v
