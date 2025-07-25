MOD=metaphlan

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        metaphlan --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        metaphlan --version

