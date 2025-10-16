MOD=metacarvel

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        python /opt/MetaCarvel/run.py --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        python /opt/MetaCarvel/run.py --help






