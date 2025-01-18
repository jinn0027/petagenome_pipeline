SINGULARITY?=/usr/local/bin/apptainer
DEF=megahit.def
SBX=megahit.sbx
SIF=megahit.sif

all : ${SBX} ${SIF}

clean :
	sudo rm -rf ${SBX} ${SIF} build-temp-*

${SBX} : ${DEF}
	sudo ${SINGULARITY} build --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	sudo ${SINGULARITY} build --fix-perms ${SIF} ${SBX}

shell : ${SBX}
	sudo ${SINGULARITY} shell --writable ${SBX} *~

help : ${SIF}
	${SINGULARITY} exec ${SIF} megahit --help

version : ${SIF}
	${SINGULARITY} exec ${SIF} megahit --version
