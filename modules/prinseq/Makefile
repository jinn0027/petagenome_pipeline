SINGULARITY?=/usr/local/bin/apptainer
DEF=prinseq.def
SBX=prinseq.sbx
SIF=prinseq.sif

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
	${SINGULARITY} exec ${SIF} prinseq-lite.pl -help

version : ${SIF}
	${SINGULARITY} exec ${SIF} prinseq-lite.pl -version

