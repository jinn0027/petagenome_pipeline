SINGULARITY?=apptainer
DEF=metaphlan2.def
SBX=metaphlan2.sbx
SIF=metaphlan2.sif

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
	${SINGULARITY} exec ${SIF} metaphlan2.py --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} metaphlan2.py --version

