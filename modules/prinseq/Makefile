MOD=prinseq

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prinseq-lite.pl -help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} prinseq-lite.pl -version

