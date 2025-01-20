SINGULARITY?=/usr/local/bin/apptainer
DEF=cdhit.def
SBX=cdhit.sbx
SIF=cdhit.sif

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
	${SINGULARITY} exec ${SIF} cd-hit --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} cd-hit --help | head -n 1

