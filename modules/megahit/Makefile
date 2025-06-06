SINGULARITY?=apptainer
DEF=megahit.def
SBX=megahit.sbx
SIF=megahit.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
	if [ -d ${SBX} ] ; then \
	    ${SINGULARITY} exec --pwd /opt --fakeroot --writable ${SBX} sh -c \
	        "rm -rf /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	fi
	rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${DEF}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${SBX}

.PHONY : shell
shell : ${SBX}
	${SINGULARITY} shell --pwd /opt --fakeroot --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} megahit --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} megahit --version
