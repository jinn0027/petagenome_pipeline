MOD=qiime2

BUILD_DEP_DIRS=../common

include ../common/common.mk

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        qiime greengenes2 non-v4-16s --help

.PHONY : lookqza
lookqza : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        -B  ${QZAFILE}:/opt/a.qza \
                        --pwd /opt \
                        ${SIF} \
                        qiime tools export \
                        --input-path /opt/a.qza \
                        --output-path /tmp/output_dir && find /tmp/output_dir -type t -print | xargs head

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --no-home \
                        --pwd /opt \
                        ${SIF} \
                        qiime --version
