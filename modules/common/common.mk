SINGULARITY?=apptainer
UNREMOVABLE_DIRS+=/usr/share/polkit-1/rules.d \
                  /etc/polkit-1/rules.d \
                  /root/go \
                  /opt

MOD?=noname
DEF=${MOD}.def
SBX=${MOD}.sbx
SIF=${MOD}.sif
OVL=${MOD}.ovl

BUILD_DEP?=
BUILD_DEP_DIRS?=

.PHONY : all
all : ${SIF}

.PHONY : clean
clean :
	@if [ -d ${SBX} ] ; then \
		if [ -d ${OVL} ] ; then \
	        for i in ${UNREMOVABLE_DIRS} ; do \
	            ${SINGULARITY} exec --no-home \
                                    --pwd /opt \
                                    --fakeroot \
                                    --overlay ${OVL} \
                                    ${SBX} \
                                    sh -c "if [ -e $${i} ] ; then chmod -R 777 $${i} ; fi"; \
	        done \
	    fi ; \
	    for i in ${UNREMOVABLE_DIRS} ; do \
	        ${SINGULARITY} exec --no-home \
                                --pwd /opt \
                                --fakeroot \
                                --writable \
                                ${SBX} \
                                sh -c "if [ -e $${i} ] ; then chmod -R 777 $${i} ; fi"; \
	    done \
	fi
	@rm -rf ${SBX} ${SIF} ${OVL} build-temp-* *~

${SIF} : ${DEF} ${BUILD_DEP}
	@for i in ${BUILD_DEP_DIRS} ; do \
	    pushd $${i}; make; popd; \
	done
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${DEF}

${SBX} : ${SIF} ${BUILD_DEP}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${SIF}

.PHONY : shell
shell : ${SBX}
	@mkdir -p ${OVL}
	${SINGULARITY} shell --no-home \
                         --pwd /opt \
                         --fakeroot \
                         --overlay ${OVL} \
                         --shell /bin/bash \
                         ${SBX}
