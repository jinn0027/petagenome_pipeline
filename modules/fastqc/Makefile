SINGULARITY?=/usr/local/bin/apptainer
DEF=fastqc.def
SBX=fastqc.sbx
SIF=fastqc.sif

all : ${SBX} ${SIF}

clean :
	sudo rm -rf ${SBX} ${SIF} build-temp-* *~

${SBX} : ${DEF}
	sudo ${SINGULARITY} build --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	sudo ${SINGULARITY} build --fix-perms ${SIF} ${SBX}

shell : ${SBX}
	sudo ${SINGULARITY} shell --writable ${SBX}

help : ${SIF}
	${SINGULARITY} exec ${SIF} fastqc --help

version : ${SIF}
	${SINGULARITY} exec ${SIF} fastqc --version
