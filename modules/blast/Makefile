SINGULARITY?=/usr/local/bin/apptainer
DEF=blast.def
SBX=blast.sbx
SIF=blast.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
	/bin/rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${DEF}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${SBX}

.PHONY : shell
shell : ${SBX}
	${SINGULARITY} shell --fakeroot --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec ${SIF} blastn -help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} blastn -version

