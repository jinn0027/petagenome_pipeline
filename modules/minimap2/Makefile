SINGULARITY?=/usr/local/bin/apptainer
DEF=minimap2.def
SBX=minimap2.sbx
SIF=minimap2.sif

.PHONY : all
all : ${SBX} ${SIF}

.PHONY : clean
clean :
	sudo rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${DEF}
	sudo ${SINGULARITY} build --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	sudo ${SINGULARITY} build --fix-perms ${SIF} ${SBX}

.PHONY : shell
shell : ${SBX}
	sudo ${SINGULARITY} shell --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec ${SIF} minimap2 2>&1

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} minimap2 2>&1 | head -n 8


