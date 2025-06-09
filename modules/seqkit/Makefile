SINGULARITY?=apptainer
DEF=seqkit.def
SBX=seqkit.sbx
SIF=seqkit.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
	if [ -d ${SBX} ] ; then \
	    ${SINGULARITY} exec --pwd /opt --fakeroot --writable ${SBX} sh -c \
	        "rm -rf /root/go /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	fi
	rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${SIF}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${SIF}

${SIF} : ${DEF}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${DEF}

.PHONY : shell
shell : ${SBX}
	${SINGULARITY} shell --pwd /opt --fakeroot --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} seqkit --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec --pwd /opt ${SIF} seqkit version

