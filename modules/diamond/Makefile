MOD=diamond

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        diamond help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        diamond version

