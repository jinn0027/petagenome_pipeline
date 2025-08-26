MOD=megan

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        blast2rma --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        blast2rma --help 2>&1 | tail -n 3

