MOD=minimap2

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} minimap2 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} minimap2 2>&1 | head -n 8


