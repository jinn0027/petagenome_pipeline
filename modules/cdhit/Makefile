MOD=cdhit

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} cd-hit --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} cd-hit --help | head -n 1

