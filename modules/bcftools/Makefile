SINGULARITY=/usr/local/bin/apptainer
DEF=bcftools.def
SBX=bcftools.sbx
SIF=bcftools.sif

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
	${SINGULARITY} exec ${SIF} bcftools --help

version : ${SIF}
	${SINGULARITY} exec ${SIF} bcftools --version

