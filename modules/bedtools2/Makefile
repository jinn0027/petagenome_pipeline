SINGULARITY?=apptainer
DEF=bedtools2.def
SBX=bedtools2.sbx
SIF=bedtools2.sif

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
	${SINGULARITY} exec ${SIF} bedtools --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} bedtools --version
