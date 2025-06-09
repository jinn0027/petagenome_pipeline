SINGULARITY?=/usr/bin/env apptainer

.PHONY : all
all : el9.sif gcc.sif java.sif golang.sif

.PHONY : sbx
sbx : el9.sbx gcc.sbx java.sbx golang.sbx

.PHONY : clean
clean :
	@for i in $$(ls | grep '.sbx') ; do \
	    ${SINGULARITY} exec --fakeroot --writable $$i sh -c \
	    "rm -rf /root/golang /usr/share/polkit-1/rules.d /etc/polkit-1/rules.d"; \
	done
	rm -rf *.sif *.sbx build-temp-* *~

### almalinux9

el9.sbx : el9.sif
	+make -f Makefile.el9 el9.sbx

el9.sif : el9.def
	+make -f Makefile.el9 el9.sif

### gcc

gcc.sbx : gcc.sif
	+make -f Makefile.gcc gcc.sbx

gcc.sif : gcc.def el9.sif
	+make -f Makefile.gcc gcc.sif

### java

java.sbx : java.sif
	+make -f Makefile.java java.sbx

java.sif : java.def el9.sif
	+make -f Makefile.java java.sif

### golang

golang.sbx : golang.sif
	+make -f Makefile.golang golang.sbx

golang.sif : golang.def el9.sif
	+make -f Makefile.golang golang.sif
