MOD=sspace

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        SSPACE_Basic.pl

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        sh -c 'ls /opt/sspace_basic/SSPACE_Basic_v*.pl | sed 's#/opt/sspace_basic/##''





