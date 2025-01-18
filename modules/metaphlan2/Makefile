SINGULARITY?=/usr/local/bin/apptainer
DEF=metaphlan2.def
SBX=metaphlan2.sbx
SIF=metaphlan2.sif

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
	${SINGULARITY} exec ${SIF} metaphlan2.py --help

version : ${SIF}
	${SINGULARITY} exec ${SIF} metaphlan2.py --version

