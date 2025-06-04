SINGULARITY?=apptainer
DEF=bwamem2.def
SBX=bwamem2.sbx
SIF=bwamem2.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
	if [ -d ${SBX} ] ; then \
	    ${SINGULARITY} exec --fakeroot --writable ${SBX} sh -c \
	        "rm -rf /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	fi
	rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${DEF}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${SBX}

.PHONY : shell
shell : ${SBX}
	${SINGULARITY} shell --fakeroot --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec ${SIF} bwa-mem2 mem 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} bwa-mem2 2>&1 | head -n 8


