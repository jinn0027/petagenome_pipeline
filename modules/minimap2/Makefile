SINGULARITY?=apptainer
DEF=minimap2.def
SBX=minimap2.sbx
SIF=minimap2.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
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
	${SINGULARITY} exec ${SIF} minimap2 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} minimap2 2>&1 | head -n 8


