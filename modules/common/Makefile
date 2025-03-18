SINGULARITY?=/usr/local/bin/apptainer

.PHONY : all
all : el9.sif gcc.sif java.sif golang.sif

.PHONY : sbx
sbx : el9.sbx gcc.sbx java.sbx golang.sbx

.PHONY : clean
clean :
	@for i in $$(ls | grep '.sbx') ; do \
	    ${SINGULARITY} exec --fakeroot --writable $$i sh -c \
	        "rm -rf /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	done
	rm -rf *.sif *.sbx build-temp-* *~

### almalinux9

el9.sbx : el9.def
	${SINGULARITY} build --fakeroot --fix-perms --sandbox el9.sbx el9.def

el9.sif : el9.sbx
	${SINGULARITY} build --fakeroot --fix-perms el9.sif el9.sbx

### gcc

gcc.sbx : gcc.def el9.sif
	${SINGULARITY} build --fakeroot --fix-perms --sandbox gcc.sbx gcc.def

gcc.sif : gcc.sbx
	${SINGULARITY} build --fakeroot --fix-perms gcc.sif gcc.sbx

### java

java.sbx : java.def el9.sif
	${SINGULARITY} build --fakeroot --fix-perms --sandbox java.sbx java.def

java.sif : java.sbx
	${SINGULARITY} build --fakeroot --fix-perms java.sif java.sbx

### golang

golang.sbx : golang.def el9.sif
	${SINGULARITY} build --fakeroot --fix-perms --sandbox golang.sbx golang.def

golang.sif : golang.sbx
	${SINGULARITY} build --fakeroot --fix-perms golang.sif golang.sbx


