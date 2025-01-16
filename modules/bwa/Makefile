SINGULARITY=/usr/local/bin/apptainer
DEF=bwa.def
SBX=bwa.sbx
SIF=bwa.sif

all : ${SBX} ${SIF}

clean :
	sudo rm -rf ${SBX} ${SIF} build-temp-*

${SBX} : ${DEF}
	sudo ${SINGULARITY} build --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	sudo ${SINGULARITY} build --fix-perms ${SIF} ${SBX}

shell : ${SBX}
	sudo ${SINGULARITY} shell --writable ${SBX}

help : ${SIF}
	${SINGULARITY} exec ${SIF} bwa 2>&1

version : ${SIF}
	${SINGULARITY} exec ${SIF} bwa 2>&1 | head -n 8


