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

.PHONY : all
all : ${SIF}

.PHONY : clean
clean :
	@if [ -d ${SBX} ] ; then \
		if [ -d ${OVL} ] ; then \
	        for i in ${UNREMOVABLE_DIRS} ; do \
	            ${SINGULARITY} exec --pwd /opt --fakeroot --overlay ${OVL} ${SBX} sh -c \
	                 "if [ -e $${i} ] ; then chmod -R 777 $${i} ; fi"; \
	        done \
	    fi ; \
	    for i in ${UNREMOVABLE_DIRS} ; do \
	        ${SINGULARITY} exec --pwd /opt --fakeroot --writable ${SBX} sh -c \
	             "if [ -e $${i} ] ; then chmod -R 777 $${i} ; fi"; \
	    done \
	fi
	@rm -rf ${SBX} ${SIF} ${OVL} build-temp-* *~

${SIF} : ${DEF} ${BUILD_DEP}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${DEF}

${SBX} : ${SIF} ${BUILD_DEP}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${SIF}

.PHONY : shell
shell : ${SBX}
	@mkdir -p ${OVL}
	${SINGULARITY} shell --pwd /opt --fakeroot --overlay ${OVL} --shell /bin/bash ${SBX}
