SINGULARITY?=apptainer
DEF=virsorter.def
SBX=virsorter.sbx
SIF=virsorter.sif

DEF_UPDATE_DB=virsorter_update_db.def
SBX_UPDATE_DB=virsorter_update_db.sbx
DB_UPDATED=virsorter-data-v2.updated.tar.gz
EXTERNAL_DB_UPDATED=../../external/${DB_UPDATED}

.PHONY : all
all : ${SBX} ${SIF} ${SBX_UPDATE_DB}

.PHONY : clean
clean :
	if [ -d ${SBX} ] ; then \
	    ${SINGULARITY} exec --fakeroot --writable ${SBX} sh -c \
	        "rm -rf /opt/VirSorter /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	fi
	if [ -d ${SBX_UPDATE_DB} ] ; then \
	    ${SINGULARITY} exec --fakeroot --writable ${SBX_UPDATE_DB} sh -c \
	        "rm -rf /opt/VirSorter /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	fi
	rm -rf ${SBX} ${SIF} ${SBX_UPDATE_DB} build-temp-* *~

${SBX} : ${DEF}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX} ${DEF}

${SIF} : ${SBX}
	${SINGULARITY} build --fakeroot --fix-perms ${SIF} ${SBX}

${SBX_UPDATE_DB} : ${DEF_UPDATE_DB}
	${SINGULARITY} build --fakeroot --fix-perms --sandbox ${SBX_UPDATE_DB} ${DEF_UPDATE_DB}

../../external/${DB_UPDATED} : ${SBX_UPDATE_DB}
	${SINGULARITY} run --fakeroot --writable ${SBX_UPDATE_DB}
	mv -f ${SBX_UPDATE_DB}/opt/VirSorter/${DB_UPDATED} ../../external/${DB_UPDATED}

.PHONY : db_updated
db_updated : ../../external/${DB_UPDATED}

.PHONY : shell
shell : ${SBX}
	${SINGULARITY} shell --fakeroot --writable ${SBX}

.PHONY : help
help : ${SIF}
	${SINGULARITY} exec ${SBX} wrapper_phage_contigs_sorter_iPlant.pl --help

.PHONY : version
version : ${SIF}
	${SINGULARITY} exec ${SIF} wrapper_phage_contigs_sorter_iPlant.pl --help | head -n 5
