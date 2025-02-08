SINGULARITY?=/usr/local/bin/apptainer
DEF=spades.def
SBX=spades.sbx
SIF=spades.sif

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
	${SINGULARITY} exec ${SIF} spades.py -h

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} spades.py -v
