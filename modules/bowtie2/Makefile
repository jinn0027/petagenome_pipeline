SINGULARITY?=/usr/local/bin/apptainer
DEF=bowtie2.def
SBX=bowtie2.sbx
SIF=bowtie2.sif

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
	${SINGULARITY} exec ${SIF} bowtie2 --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} bowtie2 --version

