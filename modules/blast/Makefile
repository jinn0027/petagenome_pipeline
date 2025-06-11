MOD=blast

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} blastn -help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} blastn -version

