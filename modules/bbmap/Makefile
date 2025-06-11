MOD=bbmap

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bbmap.sh --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} bbmap.sh --version

